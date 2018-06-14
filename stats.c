#include<stdio.h>
#include "stats.h"
#include<stdlib.h>

#define NBINS 10001
#define CFACTOR 0.0001
#define NITER_MAX 101

int run = 0;

// Design rationale: Provide general integer and double CDF for
// variables which can be updated without storing once created and generalized.

// ----- Frequency array update functions ----- //

/*
 * Updates a frequency array (bins) for a vector of double
 * values. nbins is the number of bins. X is the double vector whose
 * cdf is altered; xsize is its size. cfactor gives the conversion
 * factor to convert the double element to the appropriate (uniformly
 * sized) interval into which it must be binned.
 */
int doublelist_cdf(int* bins, int nbins, double* x, int xsize, double cfactor) {
  int i;

  for (i=0; i<xsize; i++)
    bins[(int)(x[i]/cfactor)] += 1;

  return 0;
}

// Updates an integer frequency array with a given integer as additional input.
int intlist_cdf(int* bins, int nbins, int x) {
  
  if (x < nbins && x >= 0)
    bins[x] += 1;
  else
    fprintf(stderr, "E: intlist_cdf: element %d out of binsize bounds (0-%d)\n", x, nbins);

  return 0;

}

// Used for frequency update AFTER convergence
int updatecdf_finals(struct msk_stats_t *s, double* alpha, double* beta, int numiter) {
  doublelist_cdf(s->alpha_final_cdf, NBINS, alpha, s->alpha_size, CFACTOR);
  doublelist_cdf(s->beta_final_cdf, NBINS, beta, s->beta_size, CFACTOR);
  intlist_cdf(s->numiter_cdf, NITER_MAX, numiter);
  return 0;
}

// Used for frequency update of the objective difference between global 
// and primally decomposed problem
int updatecdf_objdiff(struct msk_stats_t *s, double objdiff) {
  s->objdiff_values[run] = objdiff;
  run ++;

  // set objdiff_runs to this value, used as length of objdiff_values later on
  s->objdiff_runs = run;

  return 0;
}

// Used for interim frequency updates, when alpha's and beta's are still converging
int updatecdf_interim(struct msk_stats_t *s, double* alpha, double* beta) {
  doublelist_cdf(s->alpha_interim_cdf, NBINS, alpha, s->alpha_size, CFACTOR);
  doublelist_cdf(s->beta_interim_cdf, NBINS, beta, s->beta_size, CFACTOR);
  return 0;
}

// ------ Memory alloc and dealloc functions ------ //

// Initializes frequency arrays for quantities of interest
int alloc_cdf(struct msk_stats_t* s, int I, int J, int C, int N, int runs) {
  
  s->alpha_final_cdf = (int*)calloc(NBINS, sizeof(int));
  s->beta_final_cdf = (int*)calloc(NBINS, sizeof(int));
  s->alpha_interim_cdf = (int*)calloc(NBINS, sizeof(int));
  s->beta_interim_cdf = (int*)calloc(NBINS, sizeof(int));
  s->objdiff_values = (double*)calloc(runs, sizeof(double));
  s->numiter_cdf = (int*)calloc(NITER_MAX, sizeof(int));
  s->alpha_size = I*C*N;
  s->beta_size = I*J*C;

  if (s->alpha_final_cdf == 0 || s->beta_final_cdf == 0 || s->alpha_interim_cdf == 0 || s->beta_interim_cdf == 0 || s->objdiff_values == 0 || s->numiter_cdf == 0) {
    fprintf(stderr, "E: init_cdf: Could not initialize memory for statistics!\n");
    return 1;
  }
  
  return 0;

}

// deallocate statistics memory
int dealloc_cdf(struct msk_stats_t* s) {
  free(s->alpha_final_cdf);
  free(s->beta_final_cdf);
  free(s->alpha_interim_cdf);
  free(s->beta_interim_cdf);
  free(s->objdiff_values);
  free(s->numiter_cdf);
  return 0;
}

// ------ Functions for writing out the CDFs ------ //

int write_general_cdf(FILE* fp, int* bins, int nbins, double cfactor) {
  
  int i;
  int binsum = 0;
  for (i=0; i<nbins; i++)
    binsum += bins[i];

  int currsum = 0;
  for (i=0; i<nbins; i++) {
    currsum += bins[i];
    fprintf(fp, "%lf %lf\n", i*cfactor, (double)currsum/binsum);
  }

  return 0;
  
}

int write_objdiff_cdf(FILE* fp, double* values, int size) {

  int i;
  
  if (size < 2) {
    fprintf(stderr, "E: In writing objdiff_cdf: Too few elements in list\n");
    return 0;
  }
  
  // scan 1: determining range
  double minvalue = values[0], maxvalue = values[0];
  for (i=1; i<size; i++) {
    if (values[i] > maxvalue) maxvalue = values[i];
    if (values[i] < minvalue) minvalue = values[i];
  }

  if (minvalue <= -100.0001)
    fprintf(stderr, "Min objdiff is less than -100, something smells fishy here.\n");

  if (maxvalue - minvalue <= 0.000001)
    fprintf(stderr, "Maxvalue and minvalue are too close, maxvalue: %lf\n", maxvalue);

  if (maxvalue < 5) maxvalue = 5;
  if (minvalue > -5) minvalue = -5;

  // scan 2: uniform binning through range for CDF
  int* bins = (int*)calloc(NBINS, sizeof(int));
  for (i=0; i<size; i++) {
    int bin_index = (((values[i]-minvalue)/(maxvalue-minvalue))*(NBINS-1));
    bins[bin_index] ++;
  }
    
  // write down CDF values now
  int currsum = 0;
  for (i=0; i<NBINS; i++) {
    currsum += bins[i];
    double xval = (i*(maxvalue-minvalue)/(NBINS-1)) + minvalue;
    fprintf(fp, "%lf %lf\n", xval, (double)currsum/size);
  }

  return 0;

}

// allowing this function to stay for backward compatibility with metaparam
// Used for writing various CDFs when dealing with synthetic data.
int write_cdfs(struct msk_stats_t *s) {
  
  // Open statistics files
  FILE* fp_alpha_final = fopen("results/alpha-final.txt", "w");
  FILE* fp_beta_final = fopen("results/beta-final.txt", "w");
  FILE* fp_alpha_interim = fopen("results/alpha-interim.txt", "w");
  FILE* fp_beta_interim = fopen("results/beta-interim.txt", "w");
  FILE* fp_objdiff = fopen("results/objdiff.txt", "w");
  FILE* fp_numiter = fopen("results/numiter.txt", "w");

  // check file pointers
  fp_alpha_final = (fp_alpha_final == 0) ? stdout : fp_alpha_final;
  fp_beta_final = (fp_beta_final == 0) ? stdout : fp_beta_final;
  fp_alpha_interim = (fp_alpha_interim == 0) ? stdout : fp_alpha_interim;
  fp_beta_interim = (fp_beta_interim == 0) ? stdout : fp_beta_interim;
  fp_objdiff = (fp_objdiff == 0) ? stdout : fp_objdiff;
  fp_numiter = (fp_numiter == 0) ? stdout : fp_numiter;

  // write out CDFs
  write_general_cdf(fp_alpha_final, s->alpha_final_cdf, NBINS, CFACTOR);
  write_general_cdf(fp_beta_final, s->beta_final_cdf, NBINS, CFACTOR);
  write_general_cdf(fp_alpha_interim, s->alpha_interim_cdf, NBINS, CFACTOR);
  write_general_cdf(fp_beta_interim, s->beta_interim_cdf, NBINS, CFACTOR);
  write_objdiff_cdf(fp_objdiff, s->objdiff_values, s->objdiff_runs);
  write_general_cdf(fp_numiter, s->numiter_cdf, NITER_MAX, 1.0);

  // close open file descriptors
  if (fp_alpha_final != stdout) fclose(fp_alpha_final);
  if (fp_beta_final != stdout) fclose(fp_beta_final);
  if (fp_alpha_interim != stdout) fclose(fp_alpha_interim);
  if (fp_beta_interim != stdout) fclose(fp_beta_interim);
  if (fp_objdiff != stdout) fclose(fp_objdiff);
  if (fp_numiter != stdout) fclose(fp_numiter);

  return 0;

}

