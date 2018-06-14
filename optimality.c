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

// Do map-route iterations or one-time optimizations based on input modes.
int main(int argc, char** argv) {

  // structure definition for the problem
  struct params_t p;
  
  if (argc < 10) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file mapping-mode[1-5] routing-mode[1-3] outfile_plotstuff \n", argv[0]);
    exit(0);
  }

  // process inputs and keep a few things handy
  int ninstances;
  FILE* ftraffic;
  input_master_traffic_timeseries(&argv[1], &p, &ninstances, &ftraffic);

  // setup mapping and routing modes
  int mmode = atoi(argv[7]); // assumed (and NOT CHECKED) to be [1-5]
  int rmode = atoi(argv[8]); // assumed (and NOT CHECKED) to be [1-3]
  FILE* fp = fopen(argv[9], "w");
  if (fp == NULL) {
    printf("Could not open output file %s !\n", argv[9]);
  }
  
  printf("Mapping mode: %d Routing mode: %d\n", mmode, rmode);
  
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

  // Perform mapping initializations if
  // best response iterations are involved.
  if (mmode == 4 || mmode == 5) init_alpha(p, alpha); // round robin
  else if (mmode == 7) {
    
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
  }

  // Perform memory allocations for mapping and routing
  if (mmode == 3)  mapper = alloc_mapnode_capped_closest(p);
  else if (mmode == 6)  mapper = alloc_mapnode_95_capped_closest(p);
  else if (mmode == 4) mapper = alloc_mapnode_rtaware(p);
  else if (mmode == 5) mapper = alloc_mapnode(p);
  else if (mmode == 7) mapper = alloc_mapnode_95_rtaware(p);

  if (rmode == 1) router = alloc_dc_lat(p);
  else if (rmode == 2) router = alloc_dc_cost(p);
  else if (rmode == 3) router = alloc_dc(p);
  else if (rmode == 4) router = alloc_dc_95_lat(p);

  
  double timestamp;
  int instance;
  for (instance = 0; instance < ninstances; instance ++ ) {
    
    printf("Instance %d \n", instance);
    
    // (major) step 1. One time mapping optimizations

    if (mmode == 1) optimize_mapping_rr(p, alpha);
    
    else if (mmode == 2) optimize_mapping_closest(p, alpha);
    
    else if (mmode == 3) {
      setup_mapnode_capped_closest(p, mapper);
      optimize_mapnode_capped_closest(&env, p, mapper, alpha, &objvalue, 0);
    }

    else if (mmode == 6) {
      setup_mapnode_95_capped_closest(p, mapper);
      optimize_mapnode_95_capped_closest(&env, p, mapper, alpha, &objvalue, 0);
    }
  
    // (major) step 2. Perform one-time routing optimization
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
      
      // Perform alpha initializations here, fresh for this sample,
      // Needed if mapping mode involves 95% constraints
      // which may not be satisfied with the new traffic volumes and
      // old converged alpha values.
      
      if (mmode == 7) {	
	struct msk_problem_t* mapper2 = alloc_mapnode_95_capped_closest(p);
	setup_mapnode_95_capped_closest(p, mapper2);
	optimize_mapnode_95_capped_closest(&env, p, mapper2, alpha, &objvalue, 0);
	dealloc_mapnode_95_capped_closest(mapper2);
      }
      
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
	else if (mmode == 5) {
	  setup_mapnode(p, mapper, beta);
	  if (optimize_mapnode(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
	else if (mmode == 7) {
	  setup_mapnode_95_rtaware(p, mapper, beta);
	  if (optimize_mapnode_95_rtaware(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
      } // best response iterations

    } // end best response part.
    
    printf("\n");

    // Step 3: Record stuff into generic output file
    double objprice = compute_total_price_alphabeta(&p, alpha, beta);
    double objperfval = 0.0;

    if (rmode == 1) objperfval = objvalue / p.K;
    else if (rmode == 4) objperfval = compute_objperf_alphabeta(&p, alpha, beta);
    else if (rmode == 3) objperfval = (objvalue - objprice) / p.K;
    else {
      printf("Can't infer objperf from routing objective value; implement this function for result computation!\n");
      exit(0);
    }
    
    double wperfval = compute_perf_alphabeta(&p, alpha, beta);
    double TT = p.TT;
    double bpriceval = objprice * 1000000 / (65.8 * TT);
    double wrtt = compute_just_wrtt_alphabeta(&p, alpha, beta);

    // data format (7 fields)
    // instance_index(hr) timestamp objperf(ms) wperf(ms) wrtt(ms) tprice($/hour) bprice ($/GB) routing_objvalue_value
    fprintf(fp, "%d %lf %lf %lf %lf %lf %lf %lf\n", instance, timestamp, objperfval, wperfval, wrtt, 3600*objprice, bpriceval, objvalue);



    // Step 4: If there is a next instance, set it up
    if (instance != (ninstances-1))
      read_traffic_volumes(ftraffic, &p, &timestamp);

  } // for every instance

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

