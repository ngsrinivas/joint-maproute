#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "mosek.h"
#include "mapnode.h"
#include "global.h"
#include "dc.h"
#include "result-compute.h"
#include <math.h>
#include <string.h>

#define MAXITER 100

#define NUMRATIOS 4

// percentages for convergence
// MUST BE in descending order
double ratios[NUMRATIOS] = {0.01, 0.001, 0.0001, 0.00001};
char filenames[NUMRATIOS][250];

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

// Do map-route iterations and determine number of iterations for 
// convergence based on different measures of closeness for
// convergence.
int main(int argc, char** argv) {

  struct params_t p;
  
  if (argc < 9) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file outfiles-prefix convergence-mode\n", argv[0]);
    printf("Convergence-mode is 0 for fresh start for every instance; 1 for old set of variables\n");
    exit(0);
  }

  int ninstances;
  FILE* ftraffic;
  // read in first instance of the timeseries
  input_master_traffic_timeseries(&argv[1], &p, &ninstances, &ftraffic);

  // Populate filename buffers with names of necessary files
  // Open and keep these file descriptors ready
  int rindex;
  FILE* fps[NUMRATIOS];
  for (rindex = 0; rindex < NUMRATIOS; rindex ++) {
    sprintf(filenames[rindex], "%s-%lf", argv[7], ratios[rindex]);
    fps[rindex] = fopen(filenames[rindex], "w");
    if (fps[rindex] == NULL) {
      printf("Could not open file %s for writing!\n", filenames[rindex]);
      exit(0);
    }
  }

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
  double* X = (double*)malloc((I*J*C + I*J)*sizeof(double));
  
  struct msk_problem_t *mapper, *router, *xp;
  double objvalue, old_objvalue=0.0;
  double mapping_objvalue, global_objvalue;

  double TT;

  // Mapping and routing initializations
  mapper = alloc_mapnode(p);
  router = alloc_dc(p);
  xp = alloc_global(p);
  
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
    // keep track of the smallest ratio for which convergence has occurred.
    // to start off with, choose the largest ratio.
    int numratios_converged = 0;
    double nextratio = ratios[0];
    int numiter = 0;
    old_objvalue = 0.0;
    
    // Determine global objective value first by solving the convex version
    setup_global(p, xp);
    if (optimize_global(&env, p, xp, X, &global_objvalue, 0) > 0) {
      // optimization has failed!
      printf("Global optimization failed under instance %d\n", i);
      exit(0);
    }
    
    if (convmode == 0) 	init_alpha(p, alpha); // round robin initialization per instance if required

    // check convergence for each ratio, in an incremental fashion
    // either until all ratios have converged or the maximum number 
    // of iterations have been spent
    while (numratios_converged < NUMRATIOS && numiter < MAXITER) {
      
      if (check_convergence(objvalue, &old_objvalue, nextratio) > 0) {
	numiter ++;
	
	// setup routing optimization
	setup_dc(p, router, alpha);
	if (optimize_dc(&env, p, router, beta, &objvalue, 0) > 0) {
	  printf("Routing optimization failed under instance %d and iteration %d\n", i, numiter);
	  exit(0);
	}

	// setup mapping optimization
	setup_mapnode(p, mapper, beta);
	if (optimize_mapnode(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) {
	  printf("Mapping optimization failed under instance %d and iteration %d\n", i, numiter);
	  exit(0);
	}
	
	printf("Iteration %d\n", numiter);
	
      }
      
      else {
	// The ratio "nextratio" just converged
	// get a percentage difference from the absolute best
	double diff_from_global = compute_objdiff(mapping_objvalue, global_objvalue); 
	fprintf(fps[numratios_converged], "%d %d %lf\n", i, numiter, diff_from_global);
	printf("Ratio %lf, instance %d: %d iterations, %lf percent\n", nextratio, i, numiter, diff_from_global);
	
	// Setup for the next ratio that needs to converge
	numratios_converged ++;
	if (numratios_converged < NUMRATIOS) nextratio = ratios[numratios_converged];

      }

    } // until convergence of all ratios or numiter has exceeded threshold

    if (numiter >= MAXITER)
      printf("Instance %d needed more than %d iterations.\n", i, MAXITER);
    
    // Setup the next instance from the traffic series
    double timestamp;
    if (i != (ninstances-1))
      read_traffic_volumes(ftraffic, &p, &timestamp);
    
    printf("Finished instance %d\n", i);

  } // for each instance

  // Finally: deallocate and close stuff.
  free(alpha);
  free(beta);
  free(X);

  dealloc_mapnode(mapper);
  dealloc_dc(router);
  dealloc_global(xp);
  
  for (i=0; i<NUMRATIOS; i++)
    fclose(fps[i]);
  
  MSK_deleteenv(&env);
  
  return 0;

}
