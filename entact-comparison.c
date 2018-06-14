#include "params.h"
#include<stdio.h>
#include "inputs.h"
#include "global.h"
#include "mosek.h"
#include "helper.h"
#include "mapnode.h"
#include "result-compute.h"

// Entact comparison.
// General methodology: Run global problem on given inputs 
// to get X_ijc. These decisions are same for our scheme and
// Entact. We estimate the inaccuracy of performance prediction
// by comparing the wRTT with backbone latencies included and
// different mapping schemes, and the wRTT computed in the
// objective function of the optimization.

int main(int argc, char* argv[]) {
  
  struct params_t p;
  
  if (argc < 10) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file mapping-mode[1-3] outfile backbone-distances-file\n", argv[0]);
    exit(0);
  }
  
  // read inputs and initialize a parameter structure
  int ninstances;
  FILE* ftraffic;
  input_master_traffic_timeseries(&argv[1], &p, &ninstances, &ftraffic);
  
  // Fetch other inputs
  int mmode = atoi(argv[7]);
  FILE* outfile = fopen(argv[8], "a");
  if (outfile == NULL) {
    printf("Could not open outfile! %s \n", argv[8]);
    exit(0);
  }
  
  // initialize a mosek environment
  MSKenv_t env = NULL;
  init_mosek_env(&env);
  
  // alloc a global problem with these parameters.
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double* X = (double*)malloc((I*J*C + I*J)*sizeof(double));
  double* alpha = (double*)malloc((I*C*N + I*J)*sizeof(double));
  // backbone latencies
  double* blats = (double*)calloc(I*I, sizeof(double));
  read_backbone_lats_from_distances(argv[9], blats);
  
  struct msk_problem_t* mapper;
  struct msk_problem_t* xp = alloc_global(p);
  if (mmode == 3) mapper = alloc_mapnode_capped_closest(p);
  
  double global_objvalue, objvalue, timestamp;
  
  int i;
  for (i=0; i<ninstances; i++) {
    
    // 1. setup global problem and obtain optimal X_ijc for
    // traffic instance.
  
    setup_global(p, xp);
    if (optimize_global(&env, p, xp, X, &global_objvalue, 0)) {
      // action when optimization fails
      printf("Global optimization failed under instance %d and mapping scheme %d\n", i, mmode);
      exit(0);
    }
    
    // 2. setup mapping problem; fortunately these are one-time (ie, 
    // no convergence issues involved)
    
    if (mmode == 1) optimize_mapping_rr(p, alpha);
    
    else if (mmode == 2) optimize_mapping_closest(p, alpha);
    
    else if (mmode == 3) {
      setup_mapnode_capped_closest(p, mapper);
      if (optimize_mapnode_capped_closest(&env, p, mapper, alpha, &objvalue, 0) > 0) {
	printf("Mapping optimization failed under instance %d and mapping scheme 3\n", i);
      }
    }

    // 3. Compute and print values required for every instance of traffic series:
    double wrtt_with_backbone = compute_wrtt_with_backbone(&p, alpha, X, blats);
    double wrtt_ourscheme = compute_just_wrtt_X(&p, X);
    double totalperf = compute_perf_X(&p, X);
    // Fields: i wrtt_with_X wrtt_with_backbone_involved totalperf_X totalperf_with_backbone totalprice-per-hour(sameforboth)
    fprintf(outfile, "%d %lf %lf %lf %lf %lf\n", i, wrtt_ourscheme, wrtt_with_backbone, totalperf, totalperf - wrtt_ourscheme + wrtt_with_backbone, 3600*compute_total_price_X(&p, X));
    
    // 4. Read fresh traffic data for next iteration
    if (i != (ninstances-1))
      read_traffic_volumes(ftraffic, &p, &timestamp);
    
    printf("Finished instance %d\n", i);

  } // for each traffic-wise problem instance
  
  // Done.
  free(X);
  free(alpha);
  fclose(ftraffic);
  dealloc_global(xp);
  free(blats);
  fclose(outfile);
  if (mmode == 3) dealloc_mapnode_capped_closest(mapper);

  fprintf(stderr, "Finished all iterations\n");
  
  // delete mosek environment
  MSK_deleteenv(&env);
	 
  return 0;

}

