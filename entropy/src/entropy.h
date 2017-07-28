#include <iostream>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include "hdf5.h"

// macro definitions
#define VERSION "1.2b -- 03 March 2014"
#define MAR 'l' //character that identifies a locus in the data
#define PLOIDY 2 // number of allele copies per locus, may not work when != 2
#define NUMGENO 4 // number of genotypes at a locus, counts both
#define GEN 3 // number of genotypes for genotype likelihood format
// orderings of heterozygotes, still 2 homozygotes, but 2 heterozygotes
#define INDEXLOWERTRIPLUSDIAG (k * (k-1)) / 2 + n + k

extern gsl_rng * r;
// function declarations
using namespace std;


// struct definitions
struct mat_vec_container_int { 
  gsl_matrix_int * mat; 
  gsl_vector_int * vec; 
}; 

struct datacont{
  gsl_matrix_int * allelecount;
  gsl_matrix_int * nreads;
  double error;
  gsl_vector * errorvector;
  int glmodel;
  gsl_matrix ** genliks;
};

struct datadimcont{
  int nind;
  int nloci;
  int npop; // structure's k
  int nltpd; // length of vector to contain lower triangle plus diagonal of npop x npop matrix
  int writeaf; // bool int, write allele frequencies
  int qmatrix; // bool int, to use Q matrix rather than q vector
};

struct paramcont{
  gsl_matrix_int ** g; // genotype, by locus and ind, by allele copy
  gsl_matrix ** gSum; // running sum of genotypes
  gsl_matrix_int ** z; // ancestry, by locus and ind and allele copy
  gsl_matrix ** zSum; //running sum of ancestry, by locus and ind and pop
  gsl_matrix * p; // pop allele frequency, by locus and pop
  gsl_vector * pi; // ancestral allele freq, by locus
  gsl_vector * fst; //Fst by pop
  double alpha; // prior on population allele freq spectrum
  gsl_matrix * q; // admixture prop, by ind and pop
  gsl_matrix ** Q;
  gsl_vector * gamma; // prior on admixture proportions (for now, all equal)
  double deviance;  // deviance
};

struct auxcont{
  // convenience containers ... not parameters
  gsl_vector * npopdouble;
  gsl_vector * npopdouble2;
  gsl_vector * npopLTPDdouble; //lower triangle plus diagonal
  gsl_vector * npopLTPDdouble2; //lower triangle plus diagonal
  gsl_matrix * nindLTPDdouble;
  gsl_vector * npopsqrdouble;
  gsl_vector * nlocidouble;
  gsl_vector * ninddouble;
  gsl_vector_uint * npopuint;
  gsl_vector_uint * npopuint2;
  gsl_vector_uint * zSam;
  gsl_vector_uint * zSamsqr;
  gsl_matrix * qinit;
  double qinitscalar; // scalar for dirichlet
  double tunepi;
  double tunefst;
  double tunealpha;
  double tunegamma;
};

struct hdf5cont{
  hid_t file;
  hid_t datatype;
  hid_t stringdatatype;
  hid_t dataspaceLIG, dataspaceLIA, dataspaceLPM, dataspaceIPM,
    dataspaceIQM, dataspaceLM, dataspacePM, dataspaceM, dataspaceS;

  hid_t datasetGprob, datasetZprob, datasetP, datasetQ, datasetQmatrix, 
    datasetPi, datasetFst, datasetAlpha, datasetGamma, datasetDeviance, 
    datasetArgs;
  hid_t locusvector, indvector, popvector, scalar, indLTPDmatrix;

  /* dataset dimensions */ 
  hsize_t dimsLIG[3], dimsLIA[3], dimsLPM[3], dimsIPM[3], dimsIQM[3],  
    dimsIQ[2], dimsLM[2], dimsPM[2], dimsM[1], dimsS[1]; 
  herr_t status; /* not sure what this does, but it appears to be important?? */

  hsize_t sr3[3]; /* Start of hyperslab */
  hsize_t sr2[2];
  hsize_t sr1[1];
  hsize_t cr3[3]; /* Block count */
  hsize_t cr2[2];
  hsize_t cr1[1]; 
  hsize_t br3[3]; /* Block size */
  hsize_t br2[2]; 
  hsize_t br1[1];
};




// function definitions
void usage(char *);

int getNloci(string filename);
void getreads(string readfile, datacont * data);
int getNind(string filename);
void mcmcUpdate(datadimcont * datadim, paramcont * params, datacont * data, auxcont * auxvar);
 void updateGenotype(datadimcont * datadim, paramcont * params, 
		     datacont * data);
void initparams(datadimcont * datadim, paramcont * params, 
		datacont * data, auxcont * auxvar);
void setupParams(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void getqinit(string file, datadimcont * datadim, auxcont * auxvar);
void updateAncestry(datadimcont * datadim, paramcont * params, 
		    datacont * data, auxcont * auxvar);
void updateAncestryQ(datadimcont * datadim, paramcont * params, 
		    datacont * data, auxcont * auxvar);
void updatePopallelefreq(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updatePi(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updateFst(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updateAlpha(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updateAdmixprop(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updateAdmixpropQ(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updateGamma(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void updateGammaQ(datadimcont * datadim, paramcont * params, auxcont * auxvar);
void writeStep(hdf5cont * hdf5, paramcont * params, datadimcont * datadim,
	       auxcont * auxvar, int step, int burn, int thin);
void setuphdf5objects(hdf5cont * hdf5, datadimcont * datadim,
		      int mcmcL, int burn, int thin);
void fixdev(double * dev);
void deviance(datadimcont * datadim, paramcont * params, datacont * data,
	      auxcont * auxvar);
void updateSums(datadimcont * datadim, paramcont * params);
void writeSums(hdf5cont * hdf5, paramcont * params, datadimcont * datadim,
	       auxcont * auxvar, int mcmcL, int burn, int thin);
double fix0prob(double x);

void getdata(string readsfile, datacont * data, datadimcont * datadim, 
	     paramcont * params, auxcont * auxvar);
void updateGenotypeGl(datadimcont * datadim, paramcont * params, datacont * data);
void devianceGl(datadimcont * datadim, paramcont * params, datacont * data,
		auxcont * auxvar);
