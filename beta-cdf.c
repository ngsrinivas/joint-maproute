#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "mosek.h"
#include "mapnode.h"
#include "global.h"
#include "dc.h"
#include "result-compute.h"
#include "cdf.h"
#include <math.h>
#include <string.h>

#define MAXITER 100

#define NUMRATIOS 1

// beta CDF related parameters
// assume max 6 iterations per instance
#define BETA_NBINS 1000
#define BETA_MAXITEMS (p.total_links * p.C * ninstances * 6)

// ratio between successive objective values 
// used to detect convergence
double ratio = 0.01;

// Checking convergence using successive objective value ratios.
int check_convergence(double new_obj, double* old_obj, double ratio) {

  int not_converged = 1;
  
  if (*old_obj >= 0.000001) {
    if ((fabs(new_obj - *old_obj)/(*old_obj)) < ratio) not_converged = 0;
  }

  *old_obj = new_obj;
  
  return not_converged;
}

// Return difference of obj1 from obj2 as a percentage of obj2.
double compute_objdiff(double obj1, double obj2) {
  
  if (obj2 >= 0.000001)
    return fabs(obj1 - obj2) * 100 / obj2;
  return 100.0;
  
}

// Adds the set of beta values to the CDF for the current interim
// set of beta values.
int add_betacdf_items(struct cdf_t* beta_cdf, double* beta, struct params_t p) {
  
  int i, j, c;
  int I = p.I;
  int J = p.J;
  int C = p.C;
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p.cap[i*J + j] >= 0.000001) { // link exists
	for (c=0; c<C; c++) {
	  if (p.perf[i*J*C + j*C + c] >= 0.000001) { // client reachable
	    //cdf_additem(beta_cdf, beta[i*J*C + j*C + c]);
	    printf("%lf\n", beta[i*J*C + j*C + c]);
	  }
	}
      }
    }
  }
  
  return 0;
  
}

// Do map-route iterations and determine number of iterations for 
// convergence based on different measures of closeness for
// convergence.
int main(int argc, char** argv) {

  struct params_t p;
  struct cdf_t *beta_cdf;
  
  if (argc < 9) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file cdf-outfile convergence-mode\n", argv[0]);
    printf("Convergence-mode is 0 for fresh start for every instance; 1 for old set of variables\n");
    exit(0);
  }
  
  int ninstances;
  FILE* ftraffic;
  // read in first instance of the timeseries
  input_master_traffic_timeseries(&argv[1], &p, &ninstances, &ftraffic);
  
  // Initialize CDF structure for interim betas.
  // (not final betas)
  //beta_cdf = cdf_init_doublecdf(BETA_MAXITEMS, BETA_NBINS);
  
  // Initialize cdf file for beta
  /* FILE* cdffile = fopen(argv[7], "w"); */
  /* if (cdffile == NULL) { */
  /*   printf("Could not open cdf-outfile %s !\n", argv[7]); */
  /*   exit(0); */
  /* } */

  // convergence mode
  int convmode = atoi(argv[8]);

  // Now, returning to business as usual
  // init mosek environment
  MSKenv_t env = NULL;
  init_mosek_env(&env);

  // allocate basic structures
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double* alpha = (double*)malloc((I*C*N + I*J)*sizeof(double));
  double* beta  = (double*)malloc((I*J*C + I*J)*sizeof(double));
  
  struct msk_problem_t *mapper, *router;
  double objvalue, old_objvalue=0.0;
  double mapping_objvalue;

  double TT;

  // Mapping and routing initializations
  mapper = alloc_mapnode(p);
  router = alloc_dc(p);
  
  // If convergence mode is 1, we just initialize alpha
  // for teh first time, and the later instances just use
  // the previous values of alpha and beta to converge.
  // convergence mode 1 is something interesting to know,
  // but we'll evaluate the other mode (fresh start for every
  // instance) anyway.
  if (convmode == 1) init_alpha(p, alpha);

  int i;
  for (i=0; i<ninstances; i++) { // for each instance in traffic series
    
    // setup routing and mapping problems, and iterate until convergence
    int numratios_converged = 0;
    int numiter = 0;
    old_objvalue = 0.0;
    
    if (convmode == 0) 	init_alpha(p, alpha); // round robin initialization per instance if required

    // check convergence for each ratio, in an incremental fashion
    // either until all ratios have converged or the maximum number 
    // of iterations have been spent
    while (numratios_converged < NUMRATIOS && numiter < MAXITER) {
      
      if (check_convergence(objvalue, &old_objvalue, ratio) > 0) {
	numiter ++;
	
	// setup routing optimization
	setup_dc(p, router, alpha);
	if (optimize_dc(&env, p, router, beta, &objvalue, 0) > 0) {
	  fprintf(stderr, "Routing optimization failed under instance %d and iteration %d\n", i, numiter);
	  exit(0);
	}

	// Add items from this (interim) beta to the beta_cdf
	add_betacdf_items(beta_cdf, beta, p);

	// setup mapping optimization
	setup_mapnode(p, mapper, beta);
	if (optimize_mapnode(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) {
	  fprintf(stderr, "Mapping optimization failed under instance %d and iteration %d\n", i, numiter);
	  exit(0);
	}
	
	fprintf(stderr, "Iteration %d\n", numiter);
	
      }
      
      else {
	// It just converged.
	fprintf(stderr, "Ratio %lf, instance %d: %d iterations\n", ratio, i, numiter);
	numratios_converged ++;
	break;
      }

    } // until convergence or iterations have exceeded threshold

    if (numiter >= MAXITER)
      fprintf(stderr, "Instance %d needed more than %d iterations.\n", i, MAXITER);
    
    // Setup the next instance from the traffic series
    double timestamp;
    if (i != (ninstances-1))
      read_traffic_volumes(ftraffic, &p, &timestamp);
    
    fprintf(stderr, "Finished instance %d\n", i);

  } // for each instance
  
  // All instances done. Compute CDF for the interim betas.
  //cdf_compute_doublecdf(beta_cdf);
  //cdf_writecdf(cdffile, beta_cdf);
  
  // Finally: deallocate and close stuff.
  free(alpha);
  free(beta);

  dealloc_mapnode(mapper);
  dealloc_dc(router);
  
  //cdf_dealloc_doublecdf(beta_cdf);
  //fclose(cdffile);
  
  MSK_deleteenv(&env);
  
  return 0;

}
