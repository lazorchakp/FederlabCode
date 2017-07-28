// file: main.C

// entropy - a program not unlike structure

// Time-stamp: <Monday, 03 March 2014, 10:28 MST -- zgompert>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
//#include <getopt.h>
#include <unistd.h>
#include "hdf5.h"
#include "entropy.h"


// TODO: add storage of command line options to hdf5 object

/* To compile on OSX  */
/* h5c++ -Wall -O2 -o entropy main.C func.C -lgsl -lm  */


// compile on seismic 

// g ++ -O2 -I/usr/local/gsl-1.14/include
// -I/usr/local/atlas-3.9.28/include -L/usr/local/gsl-1.14/lib
// -L/usr/local/atlas-3.9.28/lib -Wall -o alleleEst main.C func.C
// -lgsl -lcblas -latlas -lm

using namespace std;

gsl_rng * r;  /* global state variable for random number generator */

/* ----------------- */
/* beginning of main */
/* ----------------- */

int main(int argc, char *argv[]) {

  time_t start = time(NULL);
  time_t end;
  int rng_seed = 0, ch = 0, step = 0, mcmcL= 10000;
  int thin = 1, burn = 1000;

  string readsfile = "undefined";
  string initqfile = "undefined";
  paramcont params;
  auxcont auxvar;

  hdf5cont hdf5;
  
  auxvar.qinitscalar = 1;
  auxvar.tunepi = 0.1;
  auxvar.tunefst = 0.01;
  auxvar.tunegamma = 0.2;
  auxvar.tunealpha = 0.1;
  datacont data;
  datadimcont datadim;

  datadim.npop = 2;
  datadim.writeaf = 1; // default write population allele frequencies, turn off to save disk-space and increase write speed
  datadim.qmatrix = 0; // default is to use simple P(z) = q model, rather than a matrix of two allele ancestry probabilities
  data.error = 0.0; // set to 9 for locus-specific error rates

  /* Create a new file using H5F_ACC_TRUNC access, default file
     creation properties, and default file access properties. */
  char * hdf5outfile =  (char *) "mcmcout.hdf5";

  // set up gsl random number generation 
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  srand(time(NULL)*getpid()); // VSC - Added pid to avoid equal runs when they start at the same time
  rng_seed = rand();

  // get command line arguments
  if (argc < 2) {
    usage(argv[0]);
  }
  
  while ((ch = getopt(argc, argv, "b:l:i:t:o:e:p:f:y:a:k:m:w:Q:q:s:r:")) != -1){
    switch(ch){
    case 'b':
      burn = atoi(optarg);
      break;
    case 'l':
      mcmcL = atoi(optarg);
      break;
    case 'i':
      readsfile = optarg;
      break;
    case 't':
      thin = atoi(optarg);
      break;
    case 'o':
      hdf5outfile = optarg;
      break;
    case 'e':
      data.error = atof(optarg);
      break;
    case 'p':
      auxvar.tunepi = atof(optarg);
      break;
    case 'f':
      auxvar.tunefst = atof(optarg);
      break;
    case 'y':
      auxvar.tunegamma = atof(optarg);
      break;
    case 'a':
      auxvar.tunealpha = atof(optarg);
      break;
    case 'k':
      datadim.npop = atoi(optarg);
      break;
    case 'm':
      data.glmodel = atoi(optarg);
      break;
    case 'w':
      datadim.writeaf = atoi(optarg);
      break;
    case 'Q':
      datadim.qmatrix = atoi(optarg);
      break;
    case 'q':
      initqfile = optarg;
      break;
    case 's':
      auxvar.qinitscalar = atof(optarg);
      break;
    case 'r':
      rng_seed = atoi(optarg);
      break;
    case '?':
    default:
      usage(argv[0]);
    }
  }

  gsl_rng_set(r, rng_seed); /* seed gsl_rng with output of rand, which
			       was seeded with result of time(NULL)  or from -r */
  

  datadim.nltpd = (datadim.npop * (datadim.npop - 1)) / 2 +  datadim.npop;
  cerr << "Seed: " << rng_seed << endl;
  cerr << "Reading input files" << endl;

  if (data.glmodel == 0){ // for snpcnt format
    // determine number of loci and number of alleles per locus
    datadim.nloci = getNloci(readsfile);
    datadim.nind = getNind(readsfile);
    cerr << "Number of loci: " << datadim.nloci << endl;
    cerr << "Number of individuals: " << datadim.nind << endl;
    
    setupParams(&datadim, &params, &auxvar);
    
    // allocate memory for sequence error rates
    data.errorvector = gsl_vector_calloc(datadim.nloci);

    // allocate memory for raw read data for loci and inds
    data.allelecount = gsl_matrix_int_calloc(datadim.nloci, datadim.nind);
    data.nreads = gsl_matrix_int_calloc(datadim.nloci, datadim.nind);
    auxvar.zSam = gsl_vector_uint_calloc(datadim.npop);
    auxvar.zSamsqr = gsl_vector_uint_calloc(datadim.npop * datadim.npop);

    // read in data
    getreads(readsfile, &data);
  }
  else if (data.glmodel == 1) {// genotype likelihood format
    auxvar.zSamsqr = gsl_vector_uint_calloc(datadim.npop * datadim.npop);
    getdata(readsfile, &data, &datadim, &params, &auxvar);
  }
  else {
    cerr << "Error: infile format not SNP count or genotype likelihood." << endl;
    exit(1);
  }

  hdf5.file = H5Fcreate(hdf5outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  setuphdf5objects(&hdf5, &datadim, mcmcL, burn, thin);

  /* create string will all options according command-line option or defaults */ /* <<<<< */
  // concatenate argv into one string first before write ... TODO
  hdf5.status = H5Dwrite(hdf5.datasetArgs, hdf5.stringdatatype, H5S_ALL, 
   			 H5S_ALL, H5P_DEFAULT, argv);   

  cerr << "Initializing chain" << endl;
  // read starting values for q, if available, otherwise all 1/k
  auxvar.qinit = gsl_matrix_calloc(datadim.nind, datadim.npop);
  gsl_matrix_set_all(auxvar.qinit, (double) 1/datadim.npop);
  if (initqfile.compare("undefined") != 0)
    getqinit(initqfile, &datadim, &auxvar);
  initparams(&datadim, &params, &data, &auxvar);

  // MCMC
  cerr << "Running chain ";
  for(step=0; step<mcmcL; step++){
    mcmcUpdate(&datadim, &params, &data, &auxvar);
    if( step >= burn){
      if((step-burn) % thin == 0){
	// print mcmc output
	if (data.glmodel == 0)
	  deviance(&datadim, &params, &data, &auxvar);
	else if (data.glmodel == 1)
	  devianceGl(&datadim, &params, &data, &auxvar);	  
	writeStep(&hdf5, &params, &datadim, &auxvar, step, burn, thin);
	updateSums(&datadim, &params);
      }
      if((step % 500) == 0){
	cerr << "*";
      }
    }
  }
  cerr << endl << "Writing final results" << endl;

  // write out posteriors for z, g
  writeSums(&hdf5, &params, &datadim, &auxvar, mcmcL, burn, thin);
  // close output file
  H5Fclose(hdf5.file);

  // free dynamic memory ---------------------


  // -----------------------------------------

  // prints run time
  end = time(NULL);
  cerr << "Runtime: " << (end-start)/3600 << " hr " << (end-start)%3600/60 << " min ";
  cerr << (end-start)%60 << " sec" << endl;
  return 0;
}
