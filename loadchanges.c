#include<stdio.h>
#include<stdlib.h>
#include "params.h"
#include "mosek.h"
#include "inputs.h"
#include "mapnode.h"
#include "dc.h"
#include "result-compute.h"
#include<math.h>
#include<string.h>
#include<assert.h>

#define MAXITER 100

// KNOWN BUGS
// 1. Deallocating problem setup structures (msk_problem_t's) needs to
// be done for each iteration of the scaling. This is not a serious
// problem as long as the scaling loop is commented out, as of 8th
// october '11.

// Checking convergence using objective values.
int check_convergence(double new_obj, double* old_obj) {

  int not_converged = 1;
  
  if (*old_obj >= 0.000001) {
    if ((fabs(new_obj - *old_obj)/(*old_obj)) < 0.01) not_converged = 0;
  }

  *old_obj = new_obj;
  
  return not_converged;
}

// Record original traffic volumes for scaling later on.
int record_original_volumes(struct params_t p, double* vols, double* TT) {
  
  int C = p.C;
  int N = p.N;
  
  memcpy((void*)vols, (void*)p.vol, C*N*sizeof(double));
  *TT = p.TT;
  
  return 0;
}

int set_volumes(struct params_t* p, double* vols, double TT, double scaling) {
  
  int C = p->C;
  int N = p->N;
  int c, n;
  
  for (c=0; c<C; c++) {
    for (n=0; n<N; n++) {
      p->vol[c*N + n] = scaling * vols[c*N + n];
    }
  }
  
  p->TT = TT * scaling;
  
  return 0;
  
}

// Do map-route iterations or one-time optimizations based on input modes.
int main(int argc, char** argv) {

  // attempt to get to a pareto optimal curve.
  struct params_t p;
  
  if (argc < 11) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file mapping-mode[1-5] routing-mode[1-3] outfile scale\n", argv[0]);
    exit(0);
  }

  // process inputs and keep a few things handy
  input_master_single_instance(&argv[1], &p);
  int mmode = atoi(argv[7]); // assumed (and NOT CHECKED) to be [1-6]
  int rmode = atoi(argv[8]); // assumed (and NOT CHECKED) to be [1-3]
  double scale_input = atof(argv[10]);
  FILE* fp = fopen(argv[9], "a"); // NOT CHECKED!; append mode!!
  
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

  // Record original traffic volumes somewhere; we'll be
  // using this to scale up the traffic volumes in p.
  double* vols = (double*)malloc(C*N*sizeof(double));
  double TT;
  record_original_volumes(p, vols, &TT);
  double scaling;
  
  //for (scaling = scale_input; scaling <= scale_input; scaling += 1.0) {
  scaling = scale_input;
    printf("Scaling=%lf ", scaling);
    
    // (major) step 0. Set up traffic volumes based on the scaling factor.
    set_volumes(&p, vols, TT, scaling);

    // (major) step 1. Mapping related initializations.

    if (mmode == 1) optimize_mapping_rr(p, alpha);
    
    else if (mmode == 2) optimize_mapping_closest(p, alpha);
    
    else if (mmode == 3) {
      mapper = alloc_mapnode_capped_closest(p);
      setup_mapnode_capped_closest(p, mapper);
      optimize_mapnode_capped_closest(&env, p, mapper, alpha, &objvalue, 0);
    }
    
    else if (mmode == 6) {
      mapper = alloc_mapnode_95_capped_closest(p);
      setup_mapnode_95_capped_closest(p, mapper);
      optimize_mapnode_95_capped_closest(&env, p, mapper, alpha, &objvalue, 0);
    }
  
    else if (mmode == 4) {
      mapper = alloc_mapnode_rtaware(p);
      init_alpha(p, alpha); // round robin initialization
    }
    
    else if (mmode == 7) { // routing aware 95% capped.
      // need to start with a solution for alpha where 
      // we know there's at least one feasible beta solution
      // even if that problem involves 95% caps

      // So, start by getting a feasible alpha for 95% aggregate
      // capacity. This is guaranteed to be feasible on beta even if 
      // per link 95% capacity constraints are imposed.
      mapper = alloc_mapnode_95_capped_closest(p);
      setup_mapnode_95_capped_closest(p, mapper);
      optimize_mapnode_95_capped_closest(&env, p, mapper, alpha, &objvalue, 0);
      dealloc_mapnode_95_capped_closest(mapper);
      
      // we have a feasible alpha which can give us a feasible beta
      mapper = alloc_mapnode_95_rtaware(p);
    }
    
    else if (mmode == 5) {
      mapper = alloc_mapnode(p);
      init_alpha(p, alpha); // round robin initialization
    }
    
    // (major) step 2. Routing initializations.
    
    if (rmode == 1) router = alloc_dc_lat(p);
    else if (rmode == 4) router = alloc_dc_95_lat(p);
    else if (rmode == 2) router = alloc_dc_cost(p);
    else if (rmode == 3) router = alloc_dc(p);
    
    // (major) step 3. Perform one-time routing optimization
    // or best response depending on the mapping mode used.

    if (mmode <= 3 || mmode == 6) { // one time
      if (rmode == 1) {
	setup_dc_lat(p, router, alpha);
	optimize_dc_lat(&env, p, router, beta, &objvalue, 0);
      }

      else if (rmode == 4) {
	setup_dc_95_lat(p, router, alpha);
	optimize_dc_95_lat(&env, p, router, beta, &objvalue, 0);
      }
	
      else if (rmode == 2) {
	setup_dc_cost(p, router, alpha);
	optimize_dc_cost(&env, p, router, beta, &objvalue, 0);
      }

      else if (rmode == 3) {
	setup_dc(p, router, alpha);
	optimize_dc(&env, p, router, beta, &objvalue, 0);
      }
    }
  
    else { // best response iterations.
      int numiter = 0;
      old_objvalue = 0.0;
      while (numiter < MAXITER && check_convergence(objvalue, &old_objvalue) > 0) {
      
	numiter++;
      
	if (rmode == 1) {
	  setup_dc_lat(p, router, alpha);
	  if (optimize_dc_lat(&env, p, router, beta, &objvalue, 0) > 0) break;
	}
	else if (rmode == 4) {
	  setup_dc_95_lat(p, router, alpha);
	  if (optimize_dc_95_lat(&env, p, router, beta, &objvalue, 0) > 0) break;
	}
	else if (rmode == 2) {
	  setup_dc_cost(p, router, alpha);
	  if (optimize_dc_cost(&env, p, router, beta, &objvalue, 0) > 0) break;
	}
	else if (rmode == 3) {
	  setup_dc(p, router, alpha);
	  if (optimize_dc(&env, p, router, beta, &objvalue, 0) > 0) break;
	}
      	
	printf("iteration %d ", numiter);
	
	if (mmode == 4) {
	  setup_mapnode_rtaware(p, mapper, beta);
	  if (optimize_mapnode_rtaware(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
	else if (mmode == 7) {
	  setup_mapnode_95_rtaware(p, mapper, beta);
	  if (optimize_mapnode_95_rtaware(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}

	else if (mmode == 5) {
	  setup_mapnode(p, mapper, beta);
	  if (optimize_mapnode(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
      } // best response iterations

    } // end best response part.
    
    printf("\n");
    
    // (major) step 4. record whatever is to be recorded.
    // data format: scaling wrtt wperf total_price price_per_byte routing_objective_value
    if (rmode == 1 || rmode == 3) { 
      fprintf(fp, "%lf %lf %lf %e %e %e\n", scaling, compute_just_wrtt_alphabeta(&p, alpha, beta), compute_perf_alphabeta(&p, alpha, beta), 3600*compute_total_price_alphabeta(&p, alpha, beta), 3600*compute_price_per_byte_alphabeta(&p, alpha, beta), objvalue);
      // check if the phi_ij function values in compute_perf_alphabeta are
      // computed correctly.
      double perf = compute_perf_alphabeta(&p, alpha, beta);
      double perf2 = compute_perf2_alphabeta(&p, alpha, beta);
      /* printf("perf: %lf perf2: %lf\n", perf, perf2); */
      assert(fabs(perf - perf2) <= 0.01);
    }
    
    // the following computes and prints the performance approximation
    // objective value
    // instead of the objvalue (which is just the wrtt)
    // in the case when routing doesn't involve the link penalty variables
    // at all
    else if (rmode == 4)
      fprintf(fp, "%lf %lf %lf %e %e %e\n", scaling, compute_just_wrtt_alphabeta(&p, alpha, beta), compute_perf_alphabeta(&p, alpha, beta), 3600*compute_total_price_alphabeta(&p, alpha, beta), 3600*compute_price_per_byte_alphabeta(&p, alpha, beta), p.K * compute_objperf_alphabeta(&p, alpha, beta));
    
   else if (rmode == 2) 
      fprintf(stderr, "Need to implement functionality for printing statistics for this routing mode\n");
    
    // }

  // Finally: deallocate stuff.
  fclose(fp);
  free(alpha);
  free(beta);
  
  if (mmode == 3) dealloc_mapnode_capped_closest(mapper);
  else if (mmode == 6) dealloc_mapnode_95_capped_closest(mapper);
  else if (mmode == 4) dealloc_mapnode_rtaware(mapper);
  else if (mmode == 5) dealloc_mapnode(mapper);
  else if (mmode == 7) dealloc_mapnode_95_rtaware(mapper);
  
  if (rmode == 1) dealloc_dc_lat(router);
  else if (rmode == 2) dealloc_dc_cost(router);
  else if (rmode == 3) dealloc_dc(router);
  else if (rmode == 4) dealloc_dc_95_lat(router);
  
  MSK_deleteenv(&env);

}

