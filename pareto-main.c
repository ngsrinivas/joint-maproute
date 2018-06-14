#include "params.h"
#include<stdio.h>
#include "inputs.h"
#include "global.h"
#include "mosek.h"
#include "helper.h"
#include "result-compute.h"

int run_with_K(struct params_t p, double K, double L, MSKenv_t *envp,   struct msk_problem_t* xp, double* X, double* global_objvalue) {
  p.K = K;
  p.L = L;
  setup_global(p, xp);
  if (optimize_global(envp, p, xp, X, global_objvalue, 0)) {
    // action when optimization fails
    printf("Optimization failed under K=%lf and L=%lf\n", p.K, p.L);
    exit(0);
  }

  // 7 fields computed and printed below, in order:
  // K total-perf price-per-hour price-per-hour-per-byte global-objvalue just-wrtt performance-part-of-objective
  printf("%e %lf %lf %lf %lf %lf %lf\n", p.K, compute_perf_X(&p, X), 3600*compute_total_price_X(&p, X), 3600*compute_price_per_byte_X(&p, X), *global_objvalue, compute_just_wrtt_X(&p,X), compute_objperf_X(&p,X));
  //print_link_utilizations_X(&p, X);
  
  return 0;

}

// attempt to get to a pareto optimal curve.
int main(int argc, char* argv[]) {
  
  struct params_t p;
  
  if (argc < 7) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file\n", argv[0]);
    exit(0);
  }
  
  // read inputs and initialize a parameter structure
  input_master_single_instance(&argv[1], &p);
  
  // initialize a mosek environment
  MSKenv_t env = NULL;
  init_mosek_env(&env);
  
  // alloc a global problem with these parameters.
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double* X = (double*)malloc((I*J*C + I*J)*sizeof(double));
  
  struct msk_problem_t* xp = alloc_global(p);
  
  // change parameters in p here before setting up the global problem
  // and optimizing.
  
  double global_objvalue;

  printf("Total links is %d\n", p.total_links);
  
  //run_with_K(p, 0.0, 1.0, &env, xp, X, &global_objvalue); // best cost
  run_with_K(p, 1.0, 0.0, &env, xp, X, &global_objvalue); // best performance

  // Case 2. Intermediate values of K where cost and performance tradeoff
  double iter_K;
  for (iter_K = 1e-07; iter_K <= 5e-05; iter_K *= 1.1) {
    run_with_K(p, iter_K, 1.0, &env, xp, X, &global_objvalue);
  }
    
  // Done.
  free(X);
  dealloc_global(xp);
  fprintf(stderr, "Finished all iterations\n");
  
  // delete mosek environment
  MSK_deleteenv(&env);
	 
  return 0;

}

