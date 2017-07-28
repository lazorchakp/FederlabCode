#include <iostream>
#include <sstream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <float.h>

#include <math.h>
#include "hdf5.h"

#include "entropy.h"

using namespace std;

void usage(char * name){
  fprintf(stderr,"\n%s version %s\n\n", name, VERSION); 
  fprintf(stderr, "Usage:   entropy -i infile.txt [options]\n");
  fprintf(stderr, "-i Infile with genetic data for the population\n");
  fprintf(stderr, "-l Number of MCMC steps for the analysis [default = 10000]\n");
  fprintf(stderr, "-b Discard the first n MCMC samples as a burn-in [default = 1000]\n");
  fprintf(stderr, "-t Thin MCMC samples by recording every nth value [default = 1]\n");
  fprintf(stderr, "-k Number of population clusters [default = 2]\n");
  fprintf(stderr, "-e Probability of sequence error, set to '9' for locus-specific error rates [default = 0]\n");
  fprintf(stderr, "-Q Estimate intra- and interspecific ancestry and marginal q [0 or 1, default = 0]\n");
  fprintf(stderr, "-o HDF5 format outfile with .hdf5 suffix [default = mcmcout.hdf5]\n");
  fprintf(stderr, "-m Infile is in genotype likelihood format [default = 0]\n");
  fprintf(stderr, "-w Output includes population allele frequencies [default = 1]\n");
  fprintf(stderr, "-q File with expected starting values for admixture proportions\n");
  fprintf(stderr, "-s Scalar for Dirichlet init. of q, inversly prop. to variance [default = 1]\n");
  fprintf(stderr, "-p +/- proposal for ancestral allele frequency [default = 0.1]\n");
  fprintf(stderr, "-f +/- proposal for Fst [default = 0.01]\n");
  fprintf(stderr, "-y +/- proposal for gamma [default = 0.2]\n");
  fprintf(stderr, "-a +/- proposal for alpha [default = 0.1]\n");
  fprintf(stderr, "-r INT seed for random number generator [default = clock]\n");

  exit(1);
}

/*--------- Functions to determine data dimensions and read in data -------------*/

/* determine and record the number of individuals per locus, based on readsfile */
int getNind(string filename){
  string line;
  int lociCtr = -1;
  int indCtr = 0;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      if (lociCtr > -1){
	break;
      }
      indCtr = 0;
      lociCtr++;
    }
    else { // this is a line with count data
      indCtr++;
    }
  }
  infile.close();
  return(indCtr);
}

/* determine the number of loci from reads file */
int getNloci(string filename){
  string line;
  int lociCtr = 0;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }
  // get number of loci
  while (getline(infile, line)){ // read a line of input into line
    if ( MAR == line[0]){ // MAR is the character designating a locus
      lociCtr++;
    }
  }
  infile.close();
  return(lociCtr);
}

/* determine and record the number of alleles per locus */
int getNallele(string filename, gsl_vector_int * nallele){
  string line, oneword;
  int lociCtr = -1;
  int wordCtr = 0;
  int max = 0;
  int first = 1;
  ifstream infile;

  // open file
  infile.open(filename.c_str());
  if (!infile){
    cerr << "Cannot open file " << filename << endl;
    exit(1);
  }

  // determine number of alleles per locus
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      lociCtr++;
      first = 1;

    }
    else if (first == 1){ // this is a line with count data
      first = 0;
      wordCtr = 0;
      istringstream stream(line);
      while (stream >> oneword){ // read a word at a time
	  wordCtr++;
      }
      gsl_vector_int_set(nallele, lociCtr, wordCtr);
      if (wordCtr > max){
	max = wordCtr;
      }
    }
  }
  infile.close();
  return(max);
}

/* Get sequence read data */
void getreads(string readfile, datacont * data){
  int i = -1;
  int n = 0, k = 0;
  string line, oneword;
  ifstream infile;

  infile.open(readfile.c_str());
  if (!infile){
    cerr << "Cannot open file " << readfile << endl;
    exit(1);
  }

  // read in data 
  while (getline(infile, line)){ // read a line of input into line
    if (MAR == line[0]){ // MAR is the character designating a locus
      i++; // increment locus
      if (data->error >= 1){ // locus-specific error information is
			     // included in this field, which follows
			     // locus number
	istringstream stream(line);
	stream >> oneword;
	stream >> oneword; // discard "locus" and locus number, and position
	stream >> oneword;
	stream >> oneword;
	gsl_vector_set(data->errorvector,i,atof(oneword.c_str())); 
      }
      else { // set all errors probs. equal to data->error
	gsl_vector_set(data->errorvector,i,data->error); 
      }	
      n = 0;
    }
    else { // this is a line with count data
      istringstream stream(line);
      stream >> oneword; // read a word at a time
      k = atoi(oneword.c_str());
      gsl_matrix_int_set(data->allelecount, i, n, k); // count data 
      stream >> oneword; // read a word at a time
      k=k+atoi(oneword.c_str());
      gsl_matrix_int_set(data->nreads, i, n, k);
      n++;
    }
  }
  infile.close();
}
/* memory allocation for params*/
void setupParams(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  int i;

  params->g =  new gsl_matrix_int * [datadim->nloci]; // ... could this be an unsigned char instead?
  for(i=0; i<datadim->nloci; i++){
    params->g[i] = gsl_matrix_int_calloc(PLOIDY, datadim->nind);
  }
  // ... come back later to this
  params->gSum = new gsl_matrix * [datadim->nloci];
  for(i=0; i<datadim->nloci; i++){
    params->gSum[i] = gsl_matrix_calloc(3, datadim->nind);
  }

  params->z =  new gsl_matrix_int * [datadim->nloci];
  for(i=0; i<datadim->nloci; i++){
    params->z[i] = gsl_matrix_int_calloc(PLOIDY, datadim->nind);
  }

  params->zSum = new gsl_matrix * [datadim->nloci];
  for(i=0; i<datadim->nloci; i++){
    params->zSum[i] = gsl_matrix_calloc(datadim->nind, datadim->npop);
  }

  params->p = gsl_matrix_calloc(datadim->nloci, datadim->npop);
  params->pi = gsl_vector_calloc(datadim->nloci);
  params->fst = gsl_vector_calloc(datadim->npop);
  params->q = gsl_matrix_calloc(datadim->nind, datadim->npop);
  params->Q = new gsl_matrix * [datadim->nind];
  for(i=0; i<datadim->nind; i++){
    params->Q[i] = gsl_matrix_calloc(datadim->npop, datadim->npop);
  }
  if(datadim->qmatrix == 0){
    params->gamma = gsl_vector_calloc(datadim->npop);
  }
  else{
    params->gamma = gsl_vector_calloc(datadim->nltpd);
  }
  // auxilliary variables
  auxvar->npopdouble = gsl_vector_calloc(datadim->npop);
  auxvar->npopdouble2 = gsl_vector_calloc(datadim->npop);
  auxvar->npopLTPDdouble = gsl_vector_calloc(datadim->nltpd); 
  auxvar->npopLTPDdouble2 = gsl_vector_calloc(datadim->nltpd);
  auxvar->nindLTPDdouble = gsl_matrix_calloc(datadim->nind, datadim->nltpd);
  auxvar->npopsqrdouble = gsl_vector_calloc(datadim->npop * datadim->npop);
  auxvar->nlocidouble = gsl_vector_calloc(datadim->nloci);
  auxvar->ninddouble = gsl_vector_calloc(datadim->nind);
  auxvar->npopuint = gsl_vector_uint_calloc(datadim->npop);
  auxvar->npopuint2 = gsl_vector_uint_calloc(datadim->npop);

}

/* Get data to iniialize q*/
void getqinit(string file, datadimcont * datadim, auxcont * auxvar){
  int j, k;
  string line, oneword;
  ifstream infile;

  infile.open(file.c_str());
  if (!infile){
    cerr << "Cannot open file " << file << endl;
    exit(1);
  }

  // read in q init values
  for (j=0; j<datadim->nind; j++){
    getline(infile, line); // init values for individual i
    istringstream stream(line);
    for (k=0; k<datadim->npop; k++){
      stream >> oneword;
      gsl_matrix_set(auxvar->qinit, j, k, atof(oneword.c_str()));
    }
  }
  infile.close();
}

/*---------- Functions for parameter init. and MCMC ---------------*/

/* initialize model parameters */
void initparams(datadimcont * datadim, paramcont * params, 
		datacont * data, auxcont * auxvar){
  int i, j, k, m;
  double dev;
  double sum;

  // initialize gamma
  gsl_vector_set_all(params->gamma, gsl_ran_flat(r, 0.5, 1.5));

  // initialize global ancestry, this can come from starting values
  for(i=0; i<datadim->nind; i++){
    gsl_matrix_get_row(auxvar->npopdouble2, auxvar->qinit, i);
    gsl_vector_scale(auxvar->npopdouble2, auxvar->qinitscalar);
    gsl_ran_dirichlet(r, datadim->npop, auxvar->npopdouble2->data, 
		      auxvar->npopdouble->data);
    gsl_matrix_set_row(params->q, i, auxvar->npopdouble);
    for(k=0; k<datadim->npop; k++){
      for(j=0; j<=k; j++){
	if(k == j){
	  gsl_matrix_set(params->Q[i], k, j, gsl_matrix_get(params->q, i, k) *
			 gsl_matrix_get(params->q, i, j)); 
	}
	else{
 	  gsl_matrix_set(params->Q[i], k, j, 2 * gsl_matrix_get(params->q, i, k) *
			 gsl_matrix_get(params->q, i, j)); 
	}
      }
    }
  }

  // initialize ancestral allele frequency dist.
  params->alpha = gsl_ran_flat(r, 0.5, 1.5);

  // initialize local ancestry based on global ancestry
  for(i=0; i<datadim->nloci; i++){ 
    for(j=0; j<datadim->nind; j++){ 
      for(m=0; m<PLOIDY; m++){
	gsl_matrix_get_row(auxvar->npopdouble, params->q, j);
	gsl_ran_multinomial(r, datadim->npop, 1, auxvar->npopdouble->data,
			    auxvar->npopuint->data);
	for(k=0; k<datadim->npop; k++){
	  if(gsl_vector_uint_get(auxvar->npopuint, k) == 1){
	    gsl_matrix_int_set(params->z[i], m, j, k);
	    break;
	  }
	}
      }
    }
  }

  // set all allele frequencies to 0.5, temporary
  gsl_matrix_set_all(params->p, 0.5);	  

  // initialize genotypes with gibbs update
  if (data->glmodel == 0) // snpcnt format
    updateGenotype(datadim, params, data);
  else if (data->glmodel == 1) // genotype likelihood format
    updateGenotypeGl(datadim, params, data);

  // initialize ancestral allele frequency from genotypes
  for(i=0; i<datadim->nloci; i++){ 
    sum = 0;
    for(j=0; j<datadim->nind; j++){
      sum += (gsl_matrix_int_get(params->g[i], 0, j) +
	      gsl_matrix_int_get(params->g[i], 1, j));
    }
    sum /= (datadim->nind * 2.0);
    dev = gsl_ran_beta(r, sum * datadim->nind, (1.0-sum) * datadim->nind);
    fixdev(&dev);
    gsl_vector_set(params->pi, i, dev);
  }

  // random initialization of fst
  for(i=0; i<datadim->npop; i++){ 
    gsl_vector_set(params->fst, i, gsl_ran_flat(r, 0.1, 0.4));
  }

  // re-initialize allele frequencies
  updatePopallelefreq(datadim, params, auxvar);

  // re-initialize fst
  updateFst(datadim, params, auxvar);
}

/* mcmc wrapper */
void mcmcUpdate(datadimcont * datadim, paramcont * params, datacont * data, 
		auxcont * auxvar){
  if (data->glmodel == 0) // snpcnt format
    updateGenotype(datadim, params, data); // Gibbs sampling
  else if (data->glmodel == 1) // genotype likelihood format
    updateGenotypeGl(datadim, params, data);
  if(datadim->qmatrix == 0)
    updateAncestry(datadim, params, data, auxvar);
  else if (datadim->qmatrix == 1)
    updateAncestryQ(datadim, params, data, auxvar);
  updatePopallelefreq(datadim, params, auxvar);
  updatePi(datadim, params, auxvar);
  updateFst(datadim, params, auxvar);
  updateAlpha(datadim, params, auxvar);
  if(datadim->qmatrix == 0){
    updateAdmixprop(datadim, params, auxvar);
    updateGamma(datadim, params, auxvar);
  }
  else if (datadim->qmatrix == 1){
    updateAdmixpropQ(datadim, params, auxvar);
    updateGammaQ(datadim, params, auxvar);
  }
}


/* gibbs update of genotypes */
void updateGenotype(datadimcont * datadim, paramcont * params, datacont * data){
  int i, j;

  double gProb[NUMGENO]; 
  unsigned int samGeno[NUMGENO];
  int z[2];

  for(i=0; i<datadim->nloci; i++){
    data->error = gsl_vector_get(data->errorvector, i);
    for(j=0; j<datadim->nind; j++){
      // binomial for P(X | genotype)
      gProb[0] = gsl_ran_binomial_pdf(gsl_matrix_int_get(data->allelecount, i, j),
				      data->error, 
				      gsl_matrix_int_get(data->nreads, i, j));
      // both heterozygotes
      gProb[1] = gsl_ran_binomial_pdf(gsl_matrix_int_get(data->allelecount, i, j),
				      0.5, 
				      gsl_matrix_int_get(data->nreads, i, j));
      gProb[2] = gProb[1];
      gProb[3] = gsl_ran_binomial_pdf(gsl_matrix_int_get(data->allelecount, i, j),
				      1-data->error, 
				      gsl_matrix_int_get(data->nreads, i, j));

      // get ancestry for first and second allele copies
      z[0] = gsl_matrix_int_get(params->z[i], 0, j);
      z[1] = gsl_matrix_int_get(params->z[i], 1, j);

      // binomial for P(genotype | z,p)
      //homozygote
      gProb[0] = gProb[0] * (1-gsl_matrix_get(params->p, i, z[0])) * 
	(1-gsl_matrix_get(params->p, i, z[1])) ;
      //homozygote
      gProb[3] = gProb[3] * gsl_matrix_get(params->p, i, z[0]) * 
	gsl_matrix_get(params->p, i, z[1]) ;
      //heterozygote 1, allele copy zero contains the counted allele 'p' in the pop
      gProb[1] = gProb[1] * gsl_matrix_get(params->p, i, z[0]) * 
	(1-gsl_matrix_get(params->p, i, z[1])) ;
      //heterozygote, allele copy one contains the counted allele 'p' in the pop
      gProb[2] = gProb[2] * gsl_matrix_get(params->p, i, z[1]) * 
	(1-gsl_matrix_get(params->p, i, z[0] ));

      gsl_ran_multinomial(r, NUMGENO, 1, gProb, samGeno);
      if(samGeno[0] == 1){
	gsl_matrix_int_set(params->g[i], 0, j, 0);
	gsl_matrix_int_set(params->g[i], 1, j, 0);
      }
      else if(samGeno[1] == 1){
	gsl_matrix_int_set(params->g[i], 0, j, 1);
	gsl_matrix_int_set(params->g[i], 1, j, 0);
      }
      else if(samGeno[2] == 1){
	gsl_matrix_int_set(params->g[i], 0, j, 0);
	gsl_matrix_int_set(params->g[i], 1, j, 1);
      }
      else{
	gsl_matrix_int_set(params->g[i], 0, j, 1);
	gsl_matrix_int_set(params->g[i], 1, j, 1);
      }
    }
  }
}

void updateAncestry(datadimcont * datadim, paramcont * params, 
		    datacont * data, auxcont * auxvar){

  int i, j, k, m;
  double tmpProb = 0;

  for(i=0; i<datadim->nloci; i++){
    gsl_matrix_int_set_zero(params->z[i]);
    for(j=0; j<datadim->nind; j++){
      for(m=0; m<PLOIDY; m++){
	for(k=0; k<datadim->npop; k++){
	  if(gsl_matrix_int_get(params->g[i], m, j) == 1){
	    tmpProb = gsl_matrix_get(params->p, i, k);
	  }
	  else{
	    tmpProb = 1-gsl_matrix_get(params->p, i, k);
	  }

	  gsl_vector_set(auxvar->npopdouble, k, tmpProb *
			 gsl_matrix_get(params->q, j, k));
	}
	gsl_ran_multinomial(r, datadim->npop, 1, auxvar->npopdouble->data,
			    auxvar->zSam->data);
	for(k=0; k<datadim->npop; k++){
	  if(gsl_vector_uint_get(auxvar->zSam, k) == 1){
	    gsl_matrix_int_set(params->z[i], m, j, k);
	    break;
	  }
	}
      }
    }
  }  
}

void updateAncestryQ(datadimcont * datadim, paramcont * params, 
		     datacont * data, auxcont * auxvar){

  int i, j, k, m, n, done;
  double alleleProbA = 0;

  for(i=0; i<datadim->nloci; i++){
    gsl_matrix_int_set_zero(params->z[i]);
    for(j=0; j<datadim->nind; j++){
      for(k=0; k<datadim->npop; k++){
	for(n=0; n<datadim->npop; n++){
	  alleleProbA = 1;
	  // for(m=0; m<PLOIDY; m++){
	  if(k==n){
	    for(m=0; m<PLOIDY; m++){
	      if(gsl_matrix_int_get(params->g[i], m, j) == 1){
		alleleProbA = alleleProbA * gsl_matrix_get(params->p, i, k);
	      }
	      else{
		alleleProbA = alleleProbA * (1-gsl_matrix_get(params->p, i, k));
	      }
	    }
	    gsl_vector_set(auxvar->npopsqrdouble, 
			   (k * datadim->npop) + n, alleleProbA *
			   gsl_matrix_get(params->Q[j], k, n));
	  }
	  else{
	    if(gsl_matrix_int_get(params->g[i], 0, j) == 1){
	      alleleProbA = alleleProbA * gsl_matrix_get(params->p, i, k);
	    }
	    else{
	      alleleProbA = alleleProbA * (1-gsl_matrix_get(params->p, i, k));
	    }
	    if(gsl_matrix_int_get(params->g[i], 1, j) == 1){
	      alleleProbA = alleleProbA * gsl_matrix_get(params->p, i, n);
	    }
	    else{
	      alleleProbA = alleleProbA * (1-gsl_matrix_get(params->p, i, n));
	    }

	    if(k < n){
	      gsl_vector_set(auxvar->npopsqrdouble, 
			     (k * datadim->npop) + n, alleleProbA *
			     gsl_matrix_get(params->Q[j], n, k) / 2);
	    }
	    else{
	      gsl_vector_set(auxvar->npopsqrdouble, 
			     (k * datadim->npop) + n, alleleProbA *
			     gsl_matrix_get(params->Q[j], k, n) / 2);
	    }
	  }
	}
      }
      gsl_ran_multinomial(r, datadim->npop * datadim->npop, 1, 
			  auxvar->npopsqrdouble->data,
			  auxvar->zSamsqr->data);
      for(k=0, done=0; k<datadim->npop && done<1; k++){
	for(n=0; n<datadim->npop && done<1; n++){
	  if(gsl_vector_uint_get(auxvar->zSamsqr, (k * datadim->npop + n)) == 1){
	      gsl_matrix_int_set(params->z[i], 0, j, k);
	      gsl_matrix_int_set(params->z[i], 1, j, n);
	      done = 1;
	  }
	}
      }
    }
  }  
}


/* gibbs update of pop. allele frequencies */
void updatePopallelefreq(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  int i, j, k, m;
  double cntAllele, noncntAllele, dev;

  // use auxvar->npopuint to store sum of counted allele
  // use auxvar->npopuint2 to store sum of non-counted allele

  for(i=0; i<datadim->nloci; i++){
    gsl_vector_uint_set_zero(auxvar->npopuint);
    gsl_vector_uint_set_zero(auxvar->npopuint2);
    for(j=0; j<datadim->nind; j++){
      for(m=0; m<PLOIDY; m++){ // loop across two allele copies
	k = gsl_matrix_int_get(params->z[i], m, j);
	if(gsl_matrix_int_get(params->g[i], m, j) == 1){
	  gsl_vector_uint_set(auxvar->npopuint, k, 
			      1 + gsl_vector_uint_get(auxvar->npopuint, k));
	}
	else{
	  gsl_vector_uint_set(auxvar->npopuint2, k, 
			      1 + gsl_vector_uint_get(auxvar->npopuint2, k));
	}
      }
    }
    for(k=0; k<datadim->npop; k++){
      cntAllele = gsl_vector_uint_get(auxvar->npopuint, k) +
	(-1+1/gsl_vector_get(params->fst, k))
	* gsl_vector_get(params->pi, i);
      noncntAllele = gsl_vector_uint_get(auxvar->npopuint2, k) +
	(-1+1/gsl_vector_get(params->fst, k))
	* (1 - gsl_vector_get(params->pi, i));
      // finally to the Gibbs update
      dev = gsl_ran_beta(r, cntAllele, noncntAllele);
      fixdev(&dev);
      gsl_matrix_set(params->p, i, k, dev);
    }
  }
}

void updatePi(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  // Metropolis update of pi, ancestral allele frequecies
  int i, k;
  double oldpi, newpi;
  double oldprob, newprob;
  double theta;
  double tempprob;


  for(k=0;k<datadim->npop; k++){
    theta =  (1/gsl_vector_get(params->fst, k))-1;
    gsl_vector_set(auxvar->npopdouble, k, theta);
  }

  for(i=0; i<datadim->nloci; i++){
    oldprob = 0;
    newprob = 0;
    oldpi =  gsl_vector_get(params->pi, i);
    newpi = gsl_ran_flat(r, oldpi - auxvar->tunepi,
			 oldpi+ auxvar->tunepi);
    if(newpi < 1 && newpi > 0){  // otherwise the proposal has zero probability
      for(k=0;k<datadim->npop; k++){
	theta = gsl_vector_get(auxvar->npopdouble, k);
	tempprob = gsl_ran_beta_pdf(gsl_matrix_get(params->p, i, k),
				    oldpi * theta, (1-oldpi)*theta );
	if (gsl_isinf(tempprob) == 1)
	  tempprob = DBL_MAX;
	oldprob = oldprob + log(tempprob);
	tempprob = gsl_ran_beta_pdf(gsl_matrix_get(params->p, i, k),
				    newpi * theta, (1-newpi)*theta );
	if (gsl_isinf(tempprob) == 1)
	  tempprob = DBL_MAX;
	newprob = newprob + log(tempprob);
      }

      oldprob = oldprob +log(gsl_ran_beta_pdf(oldpi, params->alpha, params->alpha));
      newprob = newprob +log(gsl_ran_beta_pdf(newpi, params->alpha, params->alpha));

      if(log(gsl_ran_flat(r, 0, 1)) < newprob - oldprob){
	gsl_vector_set(params->pi, i, newpi);
      }
    } 
  }
}

void updateFst(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  // Metropolis update of Fst
  int i, k;
  double oldfst, newfst;
  double oldprob, newprob;
  double tempprob;
  double pi;

  for(k=0;k<datadim->npop; k++){
    oldprob = 0;
    newprob = 0;
    oldfst = gsl_vector_get(params->fst, k);
    newfst = gsl_ran_flat(r, oldfst - auxvar->tunefst,
			  oldfst+ auxvar->tunefst);
    if(newfst < 1 && newfst > 0){  // otherwise the proposal has zero probability
      for(i=0; i<datadim->nloci; i++){
	pi = gsl_vector_get(params->pi, i);
	tempprob = gsl_ran_beta_pdf(gsl_matrix_get(params->p, i, k),
				    pi * (-1+1/oldfst), (1-pi)*(-1+1/oldfst) );
	if (gsl_isinf(tempprob) == 1)
	  tempprob = DBL_MAX;
	oldprob = oldprob + log(tempprob);
	tempprob = gsl_ran_beta_pdf(gsl_matrix_get(params->p, i, k),
				    pi * (-1+1/newfst), (1-pi)*(-1+1/newfst) );
	if (gsl_isinf(tempprob) == 1)
	  tempprob = DBL_MAX;
	newprob = newprob + log(tempprob);
      }
      oldprob = oldprob +log(gsl_ran_beta_pdf(oldfst, 1, 1));
      newprob = newprob +log(gsl_ran_beta_pdf(newfst, 1, 1));

      if(log(gsl_ran_flat(r, 0, 1)) < (newprob - oldprob)){
	gsl_vector_set(params->fst, k, newfst);
      }
    } 
  }
}

void updateAlpha(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  int i;
  double newalpha;
  double oldprob = 0;
  double newprob = 0;
  const double priorupperbound = 10000;
  const double priorlowerbound = DBL_MIN;

  newalpha = gsl_ran_flat(r, params->alpha - auxvar->tunealpha, 
			  params->alpha + auxvar->tunealpha );

  if(newalpha <= priorupperbound && 
     newalpha>= priorlowerbound){  // otherwise the proposal has zero probability

    for(i=0; i<datadim->nloci; i++){
      newprob = newprob + log(gsl_ran_beta_pdf(gsl_vector_get(params->pi, i), 
					       newalpha, newalpha));
      oldprob = oldprob + log(gsl_ran_beta_pdf(gsl_vector_get(params->pi, i), 
					       params->alpha, params->alpha));
    }
    if(log(gsl_ran_flat(r, 0, 1)) < newprob - oldprob){
      params->alpha = newalpha;
    }
  }
}

// Gibbs update 
void updateAdmixprop(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  int i, j, k, m;

  for(j=0; j<datadim->nind; j++){
    gsl_vector_set_zero(auxvar->npopdouble);
    for(i=0; i<datadim->nloci; i++){
      for(m=0; m<PLOIDY; m++){
	k = gsl_matrix_int_get(params->z[i], m, j);
	gsl_vector_set(auxvar->npopdouble, k, 1+ gsl_vector_get(auxvar->npopdouble, k));
      }
    }
    for(k=0; k<datadim->npop; k++){
      gsl_vector_set(auxvar->npopdouble, k, gsl_vector_get(params->gamma, k) + 
		     gsl_vector_get(auxvar->npopdouble, k));
    }
    gsl_ran_dirichlet(r,datadim->npop, auxvar->npopdouble->data, auxvar->npopdouble2->data);
    gsl_matrix_set_row(params->q, j, auxvar->npopdouble2);
  }
}


// Gibbs update 
void updateAdmixpropQ(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  int i, j, k, n, tmp;

  for(j=0; j<datadim->nind; j++){
    gsl_vector_set_zero(auxvar->npopLTPDdouble); 

    for(i=0; i<datadim->nloci; i++){
      k = gsl_matrix_int_get(params->z[i], 0, j);
      n = gsl_matrix_int_get(params->z[i], 1, j);

      if(n > k){ // make sure we're using lower triangle for storage
	tmp = n;
	n = k;
	k = tmp;
      }
      gsl_vector_set(auxvar->npopLTPDdouble, INDEXLOWERTRIPLUSDIAG, 
		     1 + gsl_vector_get(auxvar->npopLTPDdouble, INDEXLOWERTRIPLUSDIAG));
    }
  
    for(k=0; k<datadim->nltpd; k++){
      gsl_vector_set(auxvar->npopLTPDdouble, k, gsl_vector_get(params->gamma, k) + 
		     gsl_vector_get(auxvar->npopLTPDdouble, k));
    }
    gsl_ran_dirichlet(r,datadim->nltpd, auxvar->npopLTPDdouble->data, 
		      auxvar->npopLTPDdouble2->data);
    for(k=0; k<datadim->npop; k++){
      for(n=0; n<=k; n++){
	gsl_matrix_set(params->Q[j], k, n, 
		       gsl_vector_get(auxvar->npopLTPDdouble2, INDEXLOWERTRIPLUSDIAG ));
      }
    }
  }
}

void updateGamma(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  // Metropolis update of Gamma
  int j;
  double ng, oldgamma;
  double oldprob = 0;
  double newprob = 0;

  oldgamma = gsl_vector_get(params->gamma, 0);
  ng = gsl_ran_flat(r, oldgamma - auxvar->tunegamma, 
		    oldgamma + auxvar->tunegamma);

  if(ng <= 10 && ng> 0){  // otherwise the proposal has zero probability

    gsl_vector_set_all(auxvar->npopdouble, ng); // use npopdouble for the vector of newgamma
    for(j=0; j<datadim->nind; j++){
      gsl_matrix_get_row(auxvar->npopdouble2, params->q, j);
      oldprob = oldprob + gsl_ran_dirichlet_lnpdf(datadim->npop, params->gamma->data, 
						    auxvar->npopdouble2->data);
      newprob = newprob + gsl_ran_dirichlet_lnpdf(datadim->npop, auxvar->npopdouble->data,
						    auxvar->npopdouble2->data);
    }
    if(log(gsl_ran_flat(r, 0, 1)) < newprob - oldprob){
      gsl_vector_set_all(params->gamma, ng);
    }
  }
}

void updateGammaQ(datadimcont * datadim, paramcont * params, auxcont * auxvar){
  // Metropolis update of Gamma
  int j, k, n;
  double ng, oldgamma;
  double oldprob = 0;
  double newprob = 0;

  oldgamma = gsl_vector_get(params->gamma, 0);
  ng = gsl_ran_flat(r, oldgamma - auxvar->tunegamma, 
		    oldgamma + auxvar->tunegamma);

  if(ng <= 10 && ng> 0){  // otherwise the proposal has zero probability

    gsl_vector_set_all(auxvar->npopLTPDdouble, ng); // use npopLTPDdouble for the vector of newgamma
    for(j=0; j<datadim->nind; j++){
      for(k=0; k<datadim->npop; k++){
	for(n=0; n<=k; n++){
	  gsl_vector_set(auxvar->npopLTPDdouble2, INDEXLOWERTRIPLUSDIAG, 
			 gsl_matrix_get(params->Q[j], k, n));
	}
      }
      oldprob = oldprob + gsl_ran_dirichlet_lnpdf(datadim->nltpd, params->gamma->data, 
						    auxvar->npopLTPDdouble2->data);
      newprob = newprob + gsl_ran_dirichlet_lnpdf(datadim->nltpd, 
						    auxvar->npopLTPDdouble->data,
						    auxvar->npopLTPDdouble2->data);
    }
    if(log(gsl_ran_flat(r, 0, 1)) < newprob - oldprob){
      gsl_vector_set_all(params->gamma, ng);
    }
  }
}


void setuphdf5objects(hdf5cont * hdf5, datadimcont * datadim,
		      int mcmcL, int burn, int thin){

  /* HDF5, create datasets */
  /*
   * Describe the size of the array and create the data space for
   * fixed size dataset.  I reuse dataspace and datatype, because they
   * are the same for each of the parameters.
   */
  // set dimensions
  hdf5->dimsLIG[0] = datadim->nloci; 
  hdf5->dimsLIG[1] = datadim->nind; 
  hdf5->dimsLIG[2] = 3;  // 2 homozygotes and 1 heterozygote class

  hdf5->dimsLIA[0] = datadim->nloci; 
  hdf5->dimsLIA[1] = datadim->nind; 
  hdf5->dimsLIA[2] = datadim->npop ;

  hdf5->dimsLPM[0] = datadim->nloci; 
  hdf5->dimsLPM[1] = datadim->npop;
  hdf5->dimsLPM[2] = (mcmcL - burn) / thin; 

  hdf5->dimsIPM[0] = datadim->nind;
  hdf5->dimsIPM[1] = datadim->npop;
  hdf5->dimsIPM[2] = (mcmcL - burn) / thin;

  hdf5->dimsIQM[0] = datadim->nind;
  hdf5->dimsIQM[1] = datadim->nltpd;
  hdf5->dimsIQM[2] = (mcmcL - burn) / thin;

  hdf5->dimsIQ[0] = datadim->nind;
  hdf5->dimsIQ[1] = datadim->nltpd;

  hdf5->dimsLM[0] = datadim->nloci; 
  hdf5->dimsLM[1] = (mcmcL - burn) / thin; 

  hdf5->dimsPM[0] = datadim->npop; 
  hdf5->dimsPM[1] = (mcmcL - burn) / thin; 

  hdf5->dimsM[0] = (mcmcL - burn) / thin;

  hdf5->dimsS[0] = 1;

  // create dataspaces
  hdf5->dataspaceLIG = H5Screate_simple(3,  hdf5->dimsLIG, NULL); /* RANK is 3 */
  hdf5->dataspaceLIA = H5Screate_simple(3,  hdf5->dimsLIA, NULL); /* RANK is 3 */
  hdf5->dataspaceLPM = H5Screate_simple(3,  hdf5->dimsLPM, NULL); /* RANK is 3 */
  hdf5->dataspaceIPM = H5Screate_simple(3,  hdf5->dimsIPM, NULL); /* RANK is 3 */
  hdf5->dataspaceIQM = H5Screate_simple(3,  hdf5->dimsIQM, NULL); /* RANK is 3 */
  hdf5->dataspaceLM = H5Screate_simple(2,  hdf5->dimsLM, NULL); /* RANK is 2 */
  hdf5->dataspacePM = H5Screate_simple(2,  hdf5->dimsPM, NULL); /* RANK is 2 */
  hdf5->dataspaceM = H5Screate_simple(1,  hdf5->dimsM, NULL); /* RANK is 1 */
  hdf5->dataspaceS = H5Screate_simple(1, hdf5->dimsS, NULL); /* RANK is 1 */

  // define datatype by copying an existing datatype, little endian float
  hdf5->datatype = H5Tcopy(H5T_NATIVE_FLOAT); 
  // zg: changed hdf5->datatype from H5T_NATIVE_DOUBLE to reduce file size and write time
  hdf5->stringdatatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(hdf5->stringdatatype, H5T_VARIABLE);
  hdf5->status = H5Tset_order(hdf5->datatype, H5T_ORDER_LE);

  // create datasets
  hdf5->datasetGprob = H5Dcreate2(hdf5->file, "gprob", hdf5->datatype, hdf5->dataspaceLIG, 
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetZprob = H5Dcreate2(hdf5->file, "zprob", hdf5->datatype, hdf5->dataspaceLIA, 
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (datadim->writeaf == 1)
    hdf5->datasetP = H5Dcreate2(hdf5->file, "p", hdf5->datatype, hdf5->dataspaceLPM, 
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetQ = H5Dcreate2(hdf5->file, "q", hdf5->datatype, hdf5->dataspaceIPM, 
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(datadim->qmatrix == 1){
    hdf5->datasetQmatrix = H5Dcreate2(hdf5->file, "Q", hdf5->datatype, hdf5->dataspaceIQM, 
				      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  if (datadim->writeaf == 1)
    hdf5->datasetPi = H5Dcreate2(hdf5->file, "pi", hdf5->datatype, hdf5->dataspaceLM, 
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetFst = H5Dcreate2(hdf5->file, "fst", hdf5->datatype, hdf5->dataspacePM, 
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetAlpha = H5Dcreate2(hdf5->file, "alpha", hdf5->datatype, hdf5->dataspaceM, 
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetGamma = H5Dcreate2(hdf5->file, "gamma", hdf5->datatype, hdf5->dataspaceM, 
				  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetDeviance = H5Dcreate2(hdf5->file, "deviance", hdf5->datatype, hdf5->dataspaceM, 
				     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hdf5->datasetArgs = H5Dcreate2(hdf5->file, "args", hdf5->stringdatatype, hdf5->dataspaceS, 
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* need a buffer for storing a vector for each parameter
     for use in the hdf5write */
  hdf5->locusvector = H5Screate_simple(1, &hdf5->dimsLPM[0], NULL); 
  hdf5->indvector = H5Screate_simple(1, &hdf5->dimsIPM[0], NULL); 
  hdf5->popvector = H5Screate_simple(1, &hdf5->dimsPM[0], NULL); 
  hdf5->scalar = H5Screate_simple(1, &hdf5->dimsS[0], NULL);
  hdf5->indLTPDmatrix = H5Screate_simple(2, hdf5->dimsIQ, NULL);  // <<< might be a problem
}

void writeStep(hdf5cont * hdf5, paramcont * params, datadimcont * datadim,
	       auxcont * auxvar, int step, int burn, int thin){
  int i, k, n;
  double qsum;

  //  p - locus, pop, mcmc
  //  q - ind, pop, mcmc
  //  pi - locus, mcmc
  //  fst - pop, mcmc
  //  alpha - mcmc
  //  gamma - mcmc
  //  deviance - mcmc
  
  // select hyperslab for p 
  hdf5->sr3[0] = 0;
  hdf5->sr3[2] = (step - burn) / thin;
  hdf5->br3[0] = datadim->nloci;
  hdf5->br3[1] = 1;
  hdf5->br3[2] = 1;
  for(k=0; k<3; k++){
    hdf5->cr3[k] = 1;
  }

  if (datadim->writeaf == 1){
    for(k=0; k<datadim->npop; k++){
      hdf5->sr3[1] = k;
      hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceLPM, H5S_SELECT_SET, 
					 hdf5->sr3, NULL, hdf5->cr3, hdf5->br3);
      
      gsl_matrix_get_col(auxvar->nlocidouble, params->p, k);
      hdf5->status = H5Dwrite(hdf5->datasetP, H5T_NATIVE_DOUBLE, hdf5->locusvector, 
			      hdf5->dataspaceLPM, H5P_DEFAULT,  
			      auxvar->nlocidouble->data);
    }
  }


  if(datadim->qmatrix == 1){
    for(i=0; i<datadim->nind; i++){
      for(k=0; k<datadim->npop; k++){
	for(n=0; n<=k; n++){
	  gsl_matrix_set(auxvar->nindLTPDdouble, i, INDEXLOWERTRIPLUSDIAG, 
			 gsl_matrix_get(params->Q[i], k, n));
	}
      }
    }
    // select hyperslab for Q
    hdf5->br3[0] = datadim->nind;
    hdf5->br3[1] = datadim->nltpd;
    hdf5->sr3[1] = 0;
    hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceIQM, H5S_SELECT_SET, 
				       hdf5->sr3, NULL, hdf5->cr3, hdf5->br3);
      
    hdf5->status = H5Dwrite(hdf5->datasetQmatrix, H5T_NATIVE_DOUBLE, hdf5->indLTPDmatrix, 
			    hdf5->dataspaceIQM, H5P_DEFAULT,  
			    auxvar->nindLTPDdouble->data);
    // get marginal estimates for q, since we're working with Qmatrix for estimation
    for(i=0; i<datadim->nind; i++){
      for(k=0; k<datadim->npop; k++){
	for(n=0, qsum=0; n<datadim->npop; n++){
	  if(k == n) {
	    qsum = qsum + gsl_matrix_get(params->Q[i], k, n);
	  }
	  else if (k > n){
	    qsum = qsum + gsl_matrix_get(params->Q[i], k, n) / 2;
	  }
	  else {
	    qsum = qsum + gsl_matrix_get(params->Q[i], n, k) / 2;
	  }
	}
	gsl_matrix_set(params->q, i, k, qsum);
      }
    }
  }


  // select hyperslab for q
  hdf5->br3[0] = datadim->nind;
  hdf5->br3[1] = 1;
  for(k=0; k<datadim->npop; k++){
    hdf5->sr3[1] = k;
    hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceIPM, H5S_SELECT_SET, 
				       hdf5->sr3, NULL, hdf5->cr3, hdf5->br3);

    gsl_matrix_get_col(auxvar->ninddouble, params->q, k);
    hdf5->status = H5Dwrite(hdf5->datasetQ, H5T_NATIVE_DOUBLE, hdf5->indvector, 
			    hdf5->dataspaceIPM, H5P_DEFAULT,  
			    auxvar->ninddouble->data);

  }

  // select hyperslab for pi
  if (datadim->writeaf == 1){
    hdf5->sr2[0] = 0;
    hdf5->sr2[1] = (step - burn) / thin;
    hdf5->br2[0] = datadim->nloci;
    hdf5->br2[1] = 1;
    for(k=0; k<2; k++){
      hdf5->cr2[k] = 1;
    }
    hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceLM, H5S_SELECT_SET, 
				       hdf5->sr2, NULL, hdf5->cr2, hdf5->br2);
    hdf5->status = H5Dwrite(hdf5->datasetPi, H5T_NATIVE_DOUBLE, hdf5->locusvector, 
			    hdf5->dataspaceLM, H5P_DEFAULT,  
			    params->pi->data);
  }

  // select hyperslab for fst
  hdf5->sr2[0] = 0;
  hdf5->sr2[1] = (step - burn) / thin;
  hdf5->br2[0] = datadim->npop;
  hdf5->br2[1] = 1;
  for(k=0; k<2; k++){
    hdf5->cr2[k] = 1;
  }
  hdf5->status = H5Sselect_hyperslab(hdf5->dataspacePM, H5S_SELECT_SET, 
				     hdf5->sr2, NULL, hdf5->cr2, hdf5->br2);
  hdf5->status = H5Dwrite(hdf5->datasetFst, H5T_NATIVE_DOUBLE, hdf5->popvector, 
			  hdf5->dataspacePM, H5P_DEFAULT,  
			  params->fst->data);

  // select hyperslabs for alpha, gamma and deviance
  hdf5->sr1[0] = (step - burn) / thin;
  hdf5->br1[0] = 1;
  hdf5->cr1[0] = 1;

  hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceM, H5S_SELECT_SET, 
				     hdf5->sr1, NULL, hdf5->cr1, hdf5->br1);

  hdf5->status = H5Dwrite(hdf5->datasetAlpha, H5T_NATIVE_DOUBLE, hdf5->scalar, 
			  hdf5->dataspaceM, H5P_DEFAULT,  
			  &params->alpha);
  hdf5->status = H5Dwrite(hdf5->datasetGamma, H5T_NATIVE_DOUBLE, hdf5->scalar, 
			  hdf5->dataspaceM, H5P_DEFAULT,  
			  &params->gamma->data[0]);

  hdf5->status = H5Dwrite(hdf5->datasetDeviance, H5T_NATIVE_DOUBLE, hdf5->scalar, 
			  hdf5->dataspaceM, H5P_DEFAULT,  
			  &params->deviance);
}

void fixdev(double * dev){
  if(*dev == 1){
    *dev = 1-DBL_EPSILON;
  }
  else if (*dev == 0){
    *dev = DBL_EPSILON;
  }
}

void deviance(datadimcont * datadim, paramcont * params, datacont * data,
	      auxcont * auxvar){

  int i, j, m, g0, g1;
  double sumlogprob = 0.0;
  double tmpprob = 0.0;

  // sum log P( X | G)
  for(i=0; i<datadim->nloci; i++){
    data->error = gsl_vector_get(data->errorvector, i);
    for(j=0; j<datadim->nind; j++){
      g0 = gsl_matrix_int_get(params->g[i], 0, j);
      g1 = gsl_matrix_int_get(params->g[i], 1, j);
      if((g0 == 0)  & (g1==0)){
	// use fix0prob to protect against log(0), which can arise
	tmpprob = fix0prob(gsl_ran_binomial_pdf(gsl_matrix_int_get(data->allelecount, i, j),
						data->error, 
						gsl_matrix_int_get(data->nreads, i, j)));

	sumlogprob += log(tmpprob);
      }
      else if ( (g0 + g1) == 1){ // heterozygotes
	tmpprob = fix0prob(gsl_ran_binomial_pdf(gsl_matrix_int_get(data->allelecount, i, j),
						0.5, 
						gsl_matrix_int_get(data->nreads, i, j)));
	sumlogprob += log(tmpprob); 
      }
      else{ // other homozygote
	tmpprob = fix0prob(gsl_ran_binomial_pdf(gsl_matrix_int_get(data->allelecount, i, j),
						1-data->error, 
						gsl_matrix_int_get(data->nreads, i, j)));
	sumlogprob += log(tmpprob);
	
      }
    }
  }
  // sum log P( G | Z, P)
  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->nind; j++){
      for(m=0; m<PLOIDY; m++){
	if(gsl_matrix_int_get(params->g[i], m, j) == 1){
	  tmpprob = fix0prob( gsl_matrix_get(params->p, i, 
					     gsl_matrix_int_get(params->z[i], m, j)));
	  sumlogprob += log(tmpprob); 
	}
	else{
	  tmpprob = fix0prob(1-gsl_matrix_get(params->p, i, 
					      gsl_matrix_int_get(params->z[i], m, j)));
	  sumlogprob += log(tmpprob);
	}
      }
    }
  }
  params->deviance = -2 * sumlogprob;
}

void updateSums(datadimcont * datadim, paramcont * params){

  int i, j, m, k, g0, g1;

  // running sum of genotypes with each MCMC step, divide by MCMC steps to get prob
  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->nind; j++){
      g0 = gsl_matrix_int_get(params->g[i], 0, j);
      g1 = gsl_matrix_int_get(params->g[i], 1, j);
      if((g0 == 0)  & (g1==0)){
	gsl_matrix_set(params->gSum[i], 0, j, gsl_matrix_get(params->gSum[i], 0, j) + 1);
      }
      else if ( (g0 + g1) == 1){ // heterozygotes
	gsl_matrix_set(params->gSum[i], 1, j, gsl_matrix_get(params->gSum[i], 1, j) + 1);
      }
      else{ // other homozygote
	gsl_matrix_set(params->gSum[i], 2, j, gsl_matrix_get(params->gSum[i], 2, j) + 1);
      }
    }
  }

  // running sum of allele copy ancestries with each MCMC step, divide by 2*MCMC steps to get prob
  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->nind; j++){
      for(m=0; m<PLOIDY; m++){
	k = gsl_matrix_int_get(params->z[i], m, j);
	gsl_matrix_set(params->zSum[i], j, k, 1+gsl_matrix_get(params->zSum[i], j, k));
      }
    }
  }
}


void writeSums(hdf5cont * hdf5, paramcont * params, datadimcont * datadim,
	       auxcont * auxvar, int mcmcL, int burn, int thin){
  /// gSum - locus, ind, genotype
  /// zSum - locus, ind, ancestry (across both allele copies)
  int i, j, k;
  double samples;

  samples = (double) (mcmcL - burn) / thin;

  // select hyperslab for g
  for(k=0; k<3; k++){  // stride
    hdf5->cr3[k] = 1;
  }
  hdf5->br3[0] = 1; // blocksize
  hdf5->br3[1] = datadim->nind;
  hdf5->br3[2] = 1;
  for(i=0; i<datadim->nloci; i++){ //start pos
    for(j=0; j<3; j++){
      hdf5->sr3[0] = i;
      hdf5->sr3[1] = 0;
      hdf5->sr3[2] = j;
      hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceLIG, H5S_SELECT_SET, 
					 hdf5->sr3, NULL, hdf5->cr3, hdf5->br3);
      gsl_matrix_get_row(auxvar->ninddouble, params->gSum[i], j);
      gsl_vector_scale(auxvar->ninddouble, (double) 1 / samples);

      hdf5->status = H5Dwrite(hdf5->datasetGprob, H5T_NATIVE_DOUBLE, hdf5->indvector, 
			      hdf5->dataspaceLIG, H5P_DEFAULT,  
			      auxvar->ninddouble->data);

    }
  }
  // select hyperslab for z
  for(k=0; k<3; k++){  // stride
    hdf5->cr3[k] = 1;
  }
  hdf5->br3[0] = 1; // blocksize
  hdf5->br3[1] = datadim->nind;
  hdf5->br3[2] = 1;
  for(i=0; i<datadim->nloci; i++){ //start pos
    for(j=0; j<datadim->npop; j++){
      hdf5->sr3[0] = i;
      hdf5->sr3[1] = 0;
      hdf5->sr3[2] = j;
      hdf5->status = H5Sselect_hyperslab(hdf5->dataspaceLIA, H5S_SELECT_SET, 
					 hdf5->sr3, NULL, hdf5->cr3, hdf5->br3);
      gsl_matrix_get_col(auxvar->ninddouble, params->zSum[i], j);
      gsl_vector_scale(auxvar->ninddouble, (double) 1  /(2 * samples));

      hdf5->status = H5Dwrite(hdf5->datasetZprob, H5T_NATIVE_DOUBLE, hdf5->indvector, 
			      hdf5->dataspaceLIA, H5P_DEFAULT,  
			      auxvar->ninddouble->data);

    }
  }
}

double fix0prob(double x){
  if(x==0){
    x=DBL_MIN;
  }
  return(x);
}


//------- new functions for working with genotype likelihood format -----//
// function to read data, get data dimensions, and allocate data and parameter memory
void getdata(string readsfile, datacont * data, datadimcont * datadim, 
	     paramcont * params, auxcont * auxvar){
  int i, j, a;
  double genotypes[GEN];
  double gensum = 0;
  string line, element;
  ifstream infile;
  istringstream stream;

  // open genotype likelihood file
  infile.open(readsfile.c_str());
  if (!infile){
    cerr << "Cannot open file " << readsfile << endl;
    exit(1);
  }
  // read line with data dimensions
  getline(infile, line);
  stream.str(line);
  stream.clear();
  stream >> element; // number of offspring
  datadim->nind = atoi(element.c_str());  
  stream >> element; // number of loci
  datadim->nloci = atoi(element.c_str()); 

  cerr << "Number of loci: " << datadim->nloci << endl;
  cerr << "Number of individuals: " << datadim->nind << endl;

  setupParams(datadim, params, auxvar); // allocates memory, same as snpcnt version
  
  // other memory allocation, including genotype likelihoods
  data->genliks =  new gsl_matrix * [datadim->nloci];
  for(i=0; i<datadim->nloci; i++){
    data->genliks[i] = gsl_matrix_calloc(datadim->nind, GEN);
  }
  auxvar->zSam = gsl_vector_uint_calloc(datadim->npop); // needed to sample z?

  // read data
  getline(infile, line); // individual ids, these are not retained

  for(i=0; i<datadim->nloci; i++){
    getline(infile, line); // data for one locus
    stream.str(line);
    stream.clear();
    stream >> element; // locus id, this is not retained
    for(j=0; j<datadim->nind; j++){
      gensum = 0;
      for(a=0; a<GEN; a++){ // store genotype likelihoods
	stream >> element; // need to convert from phred scale: phred
	// = -10 log10(prob), prob = 10^(phred/-10)
	genotypes[a] = pow(10, (atof(element.c_str())/-10.0));
	gensum += genotypes[a];
      }
      for(a=0; a<GEN; a++){ // normalize genotype likelihoods
	genotypes[a] = genotypes[a]/gensum;
	if (genotypes[a] == 0) // set to DBL_MIN if 0
	  genotypes[a] = DBL_MIN;
	gsl_matrix_set(data->genliks[i], j, a, genotypes[a]);
      }
    }
  } 
  infile.close();
  cerr << "Finished reading genotype-likelihood format data" << endl;
}

/* gibbs update of genotypes for genotype likelihood format data*/
void updateGenotypeGl(datadimcont * datadim, paramcont * params, datacont * data){
  int i, j;

  double gProb[NUMGENO]; 
  unsigned int samGeno[NUMGENO];
  int z[2];

  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->nind; j++){
      // get gentoype likelihoods
      gProb[0] = gsl_matrix_get(data->genliks[i], j, 0); // 1st homozygote
      gProb[1] = gsl_matrix_get(data->genliks[i], j, 1);
      gProb[2] = gProb[1];
      gProb[3] = gsl_matrix_get(data->genliks[i], j, 2); // 2nd homozygote
      
      // get ancestry for first and second allele copies
      z[0] = gsl_matrix_int_get(params->z[i], 0, j);
      z[1] = gsl_matrix_int_get(params->z[i], 1, j);

      // binomial for P(genotype | z,p)
      //homozygote
      gProb[0] = gProb[0] * (1-gsl_matrix_get(params->p, i, z[0])) * 
	(1-gsl_matrix_get(params->p, i, z[1])) ;
      //homozygote
      gProb[3] = gProb[3] * gsl_matrix_get(params->p, i, z[0]) * 
	gsl_matrix_get(params->p, i, z[1]) ;
      //heterozygote 1, allele copy zero contains the counted allele 'p' in the pop
      gProb[1] = gProb[1] * gsl_matrix_get(params->p, i, z[0]) * 
	(1-gsl_matrix_get(params->p, i, z[1])) ;
      //heterozygote, allele copy one contains the counted allele 'p' in the pop
      gProb[2] = gProb[2] * gsl_matrix_get(params->p, i, z[1]) * 
	(1-gsl_matrix_get(params->p, i, z[0] ));

      gsl_ran_multinomial(r, NUMGENO, 1, gProb, samGeno);
      if(samGeno[0] == 1){
	gsl_matrix_int_set(params->g[i], 0, j, 0);
	gsl_matrix_int_set(params->g[i], 1, j, 0);
      }
      else if(samGeno[1] == 1){
	gsl_matrix_int_set(params->g[i], 0, j, 1);
	gsl_matrix_int_set(params->g[i], 1, j, 0);
      }
      else if(samGeno[2] == 1){
	gsl_matrix_int_set(params->g[i], 0, j, 0);
	gsl_matrix_int_set(params->g[i], 1, j, 1);
      }
      else{
	gsl_matrix_int_set(params->g[i], 0, j, 1);
	gsl_matrix_int_set(params->g[i], 1, j, 1);
      }
    }
  }
}

// calculte deviance with genotype likelihood format data
void devianceGl(datadimcont * datadim, paramcont * params, datacont * data,
		auxcont * auxvar){

  int i, j, m, g0, g1;
  double sumlogprob = 0.0;
  double tmpprob = 0.0;

  // sum log P( X | G)
  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->nind; j++){
      g0 = gsl_matrix_int_get(params->g[i], 0, j);
      g1 = gsl_matrix_int_get(params->g[i], 1, j);
      if((g0 == 0)  & (g1==0)){
	// use fix0prob to protect against log(0), which can arise
	tmpprob = fix0prob(gsl_matrix_get(data->genliks[i], j, 0));
	sumlogprob += log(tmpprob);
      }
      else if ( (g0 + g1) == 1){ // heterozygotes
	tmpprob = fix0prob(gsl_matrix_get(data->genliks[i], j, 1));
	sumlogprob += log(tmpprob); 
      }
      else{ // other homozygote
	tmpprob = fix0prob(gsl_matrix_get(data->genliks[i], j, 2));
	sumlogprob += log(tmpprob);
	
      }
    }
  }
  // sum log P( G | Z, P)
  for(i=0; i<datadim->nloci; i++){
    for(j=0; j<datadim->nind; j++){
      for(m=0; m<PLOIDY; m++){
	if(gsl_matrix_int_get(params->g[i], m, j) == 1){
	  tmpprob = fix0prob( gsl_matrix_get(params->p, i, 
					     gsl_matrix_int_get(params->z[i], m, j)));
	  sumlogprob += log(tmpprob); 
	}
	else{
	  tmpprob = fix0prob(1-gsl_matrix_get(params->p, i, 
					      gsl_matrix_int_get(params->z[i], m, j)));
	  sumlogprob += log(tmpprob);
	}
      }
    }
  }
  params->deviance = -2 * sumlogprob;
}
