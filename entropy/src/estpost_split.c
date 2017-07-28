/*
  Peter Lazorchak
  6/26/17
  lazorchakp@gmail.com

  Modified from source code 
*/

/* read hdf5 and write data summary */
/* Time-stamp: <Tuesday, 24 September 2013, 14:17 MDT -- zgompert> */

/* Compilation for Linux */
/* h5cc -Wall -O3 -o estpost.entropy estpost_h5_entropy.c -lgsl -lgslcblas */
/* Compilation for OSX */
/* h5cc -Wall -O3 -o estpost.entropy estpost_h5_entropy.c -lgsl -lm */


/*
  This program prints to stdout the average log probability of all chains for
  each input file. This is designed to work with the perl script clump_gen.pl.
  Modifications made to this program will affect clump_gen.pl - be careful.
*/

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_histogram.h>

#include "hdf5.h"
#include <getopt.h>

/* header, function declarations */
#define MAXFILEN 20 /* maximum number of infiles */
void usage(char * name);

void estpost(hid_t * file, int burn, int nchains);

void processHDF5ind(hid_t dataspace, hid_t * dataset, gsl_vector * onesample,
    int nind, int npop, int nsamples, int burn, int nchains);

void processHDF5prob(hid_t dataspace, hid_t * dataset, gsl_vector * onesample,
    int nsamples, int burn, int nchains);

FILE * outfpInds;

int main (int argc, char **argv) {
  int ch = 0;
  int burn = 0; /* discard the first burn samples as a burn-in */
  char * infile = "undefined"; /* filename */
  char * outfile = "postout_split.csv";

  int nchains = 0; /* number of mcmc chains = number of infiles */

  /* variables for getopt_long */
  static struct option long_options[] = {
    {"version", no_argument, 0, 'v'},
    {0, 0, 0, 0} 
  };
  int option_index = 0;

  hid_t file[MAXFILEN]; /* file handle */
  
  /*  get command line arguments */
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt_long(argc, argv, "o:b:", long_options, &option_index)) != -1){
    switch(ch){
    case 'o':
      outfile = optarg;
      break;
    case 'b':
      burn = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }

  /* open the h5 files, read only */
  while (optind < argc){
    infile = argv[optind];
    fprintf(stderr, "file = %s\n", infile);
    file[nchains] = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    ++nchains;
    ++optind;
  }

  // added this, names the output file
  outfpInds = fopen(outfile, "w");
  if (!outfpInds){
    fprintf(stderr, "Can't open %s for writing!\n", outfile);
    exit(1);
  }

  fprintf(outfpInds, "ind_pop"); // added
  int i;
  for (i = 1; i <= nchains; ++i) {
    fprintf(outfpInds, ",mean%d", i);
  }
  fprintf(outfpInds, "\n");
  
  /* main function */
  estpost(file, burn, nchains);

  fclose(outfpInds); // added this
  int c;
  for(c = 0; c < nchains; ++c) {
    H5Fclose(file[c]);
  }
  return 0;
}

/* ---------------- Functions ------------------ */

/* Prints usage */
void usage(char * name){
  fprintf(stderr,"\n%s\n\n", name); 
  fprintf(stderr, "Usage:   estpost [options] infile1.hdf5 infile2.hdf5\n");
  fprintf(stderr, "-o     Outfile [default = postout_split.csv]\n");
  fprintf(stderr, "-b     Number of additinal MCMC samples to discard for burn-in [default = 0]\n\n");
  exit(1);
}

void estpost(hid_t * file, int burn, int nchains) {
  gsl_vector * onesample;
  hid_t dataset[MAXFILEN], dataspace;
  hsize_t dims[3]; /* note 3 dimensions is maximum in this application */
  /* dimensions for mcmc samples */
  int nind = 0, npop = 0, nsamples = 0;
  int chain;

  /* we already have checked that param should exist in input hdf5 file*/
  for (chain=0; chain<nchains; ++chain){
    dataset[chain] = H5Dopen2(file[chain], "q", H5P_DEFAULT);
  }
  dataspace = H5Dget_space(dataset[0]); // we assume that all chains have the same dimensions
  H5Sget_simple_extent_ndims(dataspace);
  H5Sget_simple_extent_dims(dataspace, dims, NULL);

  /* interpret dimensions based on parameter, then calculate (if
     necessary) and print desired quantities */
  nind = dims[0];
  npop = dims[1];
  nsamples = dims[2];
  if (burn >= nsamples){
    fprintf(stderr, "Burnin exceeds number of samples\n");
    exit(1);
  }
  onesample = gsl_vector_calloc(nsamples - burn);
  fprintf(stderr, "parameter dimensions for q: ind = %d, populations = %d, samples = %d, chains = %d\n", 
   nind, npop, nsamples, nchains);

  processHDF5ind(dataspace, dataset, onesample, nind, npop, nsamples, burn, nchains);
  
  // close the dataset for q and open the data set for deviance
  for (chain=0; chain<nchains; ++chain){
    H5Dclose(dataset[chain]);
    dataset[chain] = H5Dopen2(file[chain], "deviance", H5P_DEFAULT);
  }

  dataspace = H5Dget_space(dataset[0]);
  H5Sget_simple_extent_ndims(dataspace);
  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  nsamples = dims[0];

  processHDF5prob(dataspace, dataset, onesample, nsamples, burn, nchains);
}


/* calculate summaries for parameters indexed by ind and
   population */
void processHDF5ind(hid_t dataspace, hid_t * dataset,	gsl_vector * onesample,
  int nind, int npop, int nsamples, int burn, int nchains){
  int i, j, c;
  hid_t mvector;
  hsize_t start[3];  /* Start of hyperslab */
  hsize_t count[3] = {1,1,1};  /* Block count */
  hsize_t block[3] = {1, 1, nsamples - burn};

  hsize_t samdim[1];

  /* create vector for buffer */
  samdim[0] = nsamples - burn;
  mvector = H5Screate_simple(1, samdim, NULL);

  start[2] = burn;
  /* loop through population, then locus */
  for (j = 0; j < npop; ++j) {
    for (i = 0; i < nind; ++i) {
      fprintf(outfpInds, "q_ind_%d_pop_%d", i, j);
      start[0] = i; start[1] = j; /* start[2] = burn; */
      for (c = 0; c < nchains; ++c){
	      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, block);
	      H5Dread(dataset[c], H5T_NATIVE_DOUBLE, mvector, dataspace, H5P_DEFAULT, 
		      onesample->data);
        fprintf(outfpInds, ",%.6f", gsl_stats_mean(onesample->data, 1, nsamples - burn));
      }
      fprintf(outfpInds, "\n");
    }
  } 
  H5Sclose(mvector);
}

void processHDF5prob(hid_t dataspace, hid_t * dataset, gsl_vector * onesample,
    int nsamples, int burn, int nchains) {
  int c;
  hid_t mvector;
  hsize_t start[1];  /* Start of hyperslab */
  hsize_t count[1] = {1};  /* Block count */
  hsize_t block[1];
  hsize_t samdim[1];

  /* create vector for buffer */
  samdim[0] = nsamples - burn;
  mvector = H5Screate_simple(1, samdim, NULL);
  /* loop through parameters */

  start[0] = burn;
  block[0] = (nsamples - burn);
  for (c=0; c<nchains; ++c){
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, block);
    H5Dread(dataset[c], H5T_NATIVE_DOUBLE, mvector, dataspace, H5P_DEFAULT, onesample->data);
    // each column contains the deviance of every recorded iteration
    // average them, convert to log prob, and print to stdout
    printf("%.2f\n", -0.5 * gsl_stats_mean(onesample->data, 1, nsamples - burn));
  }
  H5Sclose(mvector);
}