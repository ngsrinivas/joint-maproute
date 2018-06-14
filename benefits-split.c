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
#define COSTTOL 0.00001
#define LATSTOL 0.00001

// mapping and routing modes of importance for comparison
// number of other modes to compare
#define NUMCOMPARE 7
int int_mmodes[NUMCOMPARE] = {5, 1, 2, 3, 4, 6, 6};
int int_rmodes[NUMCOMPARE] = {3, 1, 1, 1, 1, 1, 3};
int our_mmode     = 5;
int our_rmode     = 3;

// Checking convergence using objective values.
int check_convergence(double new_obj, double* old_obj) {

  int not_converged = 1;
  
  if (*old_obj >= 0.000001) {
    if ((fabs(new_obj - *old_obj)/(*old_obj)) < 0.01) not_converged = 0;
  }

  *old_obj = new_obj;
  
  return not_converged;
}

int init_counts_diffs(double* count_cli, double* count_vol, double* difflats_cli, double* diffcost_cli, double* difflats_vol, double* diffcost_vol) {
  
  int i;
  for (i=0; i<5; i++) {
    count_cli[i] = count_vol[i] = difflats_cli[i] = diffcost_cli[i] = difflats_vol[i] = diffcost_vol[i] = 0.0;
  }
  
  return 0;
  
}

int print_diff_stats(int mmode, int rmode, double* count_cli, double* count_vol, double* diffcost_cli, double* diffcost_vol, double* difflats_cli, double* difflats_vol) {
  
  int i;
  printf("Mapping mode: %d Routing mode: %d\n", mmode, rmode);
  for (i=0; i<5; i++) {
    printf("%d %lf %lf %lf %lf %lf %lf\n", i, count_cli[i], count_vol[i], diffcost_cli[i], diffcost_vol[i], difflats_cli[i], difflats_vol[i]);
  }
  
}

int print_traffic_volumes(struct params_t p) {
  
  int c, n;
  int C = p.C;
  int N = p.N;
  
  for (c=0; c<C; c++) {
    for (n=0; n<N; n++) {
      printf("c %d n %d vol %lf\n", c, n, p.vol[c*N + n]);
    }
  }
  
  return 0;
  
}

// Do map-route iterations or one-time optimizations based on input modes.
int main(int argc, char** argv) {

  // structure definition for the problem
  struct params_t p;
  
  if (argc < 8) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file instance_number \n", argv[0]);
    printf("instance_number is zero indexed in the traffic file.\n");
    exit(0);
  }

  // process inputs and keep a few things handy
  int ninstances;
  FILE* ftraffic;
  input_master_traffic_timeseries(&argv[1], &p, &ninstances, &ftraffic);

  // setup mapping and routing modes
  int mmode;
  int rmode;
  int reqd_instance = atoi(argv[7]);
  
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
  double* int_alpha = (double*)malloc((I*C*N + I*J)*sizeof(double));
  double* int_beta  = (double*)malloc((I*J*C + I*J)*sizeof(double));
  
  struct msk_problem_t *mapper, *router;
  double objvalue, old_objvalue=0.0;
  double mapping_objvalue;

  double timestamp;
  int instance;

  // get the right problem instance to work on into 'p'
  for (instance = 0; instance < ninstances; instance ++ ) {
    
    printf("Instance %d \n", instance);
    if (instance == reqd_instance) break;
    
    // Step 4: If there is a next instance, set it up
    if (instance != (ninstances-1))
      read_traffic_volumes(ftraffic, &p, &timestamp);

  } // for every instance
  
  // diagnostic to check if traffic volumes were read correctly
  //print_traffic_volumes(p);

  int i, c;
  double *actual_alpha, *actual_beta;
  // storage for average latency and costs for the `interesting` scheme
  double* avglats = (double*)malloc(p.C * sizeof(double));
  double* avgcost = (double*)malloc(p.C * sizeof(double));

  // define a few statistics variables here:
  // average values for the current client with current scheme
  double avglats_c, avgcost_c; 
  // counts and volume weighted percentages for various entries in the benefits table
  double count_cli[5], count_vol[5], difflats_cli[5], diffcost_cli[5], difflats_vol[5], diffcost_vol[5];
  
  int n;

  for (i=0; i<NUMCOMPARE; i++) {
    
    printf("On scheme %d\n", i);
    
    mmode = int_mmodes[i];
    rmode = int_rmodes[i];
    
    // step (-1). Manage alpha and beta buffers
    // in order to retain rest of code unmodified
    // with 'alpha' calls in there.
    
    if (i == 0) {
      actual_alpha = alpha;
      actual_beta = beta;
      alpha = int_alpha;
      beta = int_beta;
    }
    
    else if (i == 1) {
      alpha = actual_alpha;
      beta = actual_beta;
    }
    
    // major step 0. Perform memory allocations for mapping and routing
    
    if (mmode == 3) mapper = alloc_mapnode_capped_closest(p);
    else if (mmode == 6) mapper = alloc_mapnode_95_capped_closest(p);
    else if (mmode == 4) mapper = alloc_mapnode_rtaware(p);
    else if (mmode == 5) mapper = alloc_mapnode(p);
    else if (mmode == 7) mapper = alloc_mapnode_95_rtaware(p);

    if (rmode == 1) router = alloc_dc_lat(p);
    else if (rmode == 4) router = alloc_dc_95_lat(p);
    else if (rmode == 2) router = alloc_dc_cost(p);
    else router = alloc_dc(p);

    // Perform mapping round robin initialization if
    // best response iterations are involved.
    // including an assertion here when it is re-run so that 
    // mapping mode 7 initializations for feasibility can be taken care
    // of later on. Note that these have to go inside the loop
    // ie, using previous alpha's doesn't guarantee feasibility in
    // any meaningful way.
    fprintf(stderr, "Need to fix alpha initializations first.\n");
    exit(0);
    if (mmode == 4 || mmode == 5 || mmode == 7) init_alpha(p, alpha);

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
      while (numiter < MAXITER && check_convergence(objvalue, &old_objvalue) > 0) {
      
	numiter++;
      
	if (rmode == 1) {
	  setup_dc_lat(p, router, alpha);
	  if (optimize_dc_lat(&env, p, router, beta, &objvalue, 0) > 0) break;
	}
	if (rmode == 4) {
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
      	
	printf("iteration %d\n", numiter);
	
	if (mmode == 4) {
	  setup_mapnode_rtaware(p, mapper, beta);
	  if (optimize_mapnode_rtaware(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
	else if (mmode == 7) {
	  setup_mapnode_rtaware(p, mapper, beta);
	  if (optimize_mapnode_rtaware(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
	else if (mmode == 5) {
	  setup_mapnode(p, mapper, beta);
	  if (optimize_mapnode(&env, p, mapper, alpha, &mapping_objvalue, 0) > 0) break;
	}
      } // best response iterations

    } // end best response part.
    
    printf("\n");    
    
    // major step 3. Deallocate mapping and routing related buffers
    if (mmode == 3) dealloc_mapnode_capped_closest(mapper);
    else if (mmode == 6) dealloc_mapnode_95_capped_closest(mapper);
    else if (mmode == 4) dealloc_mapnode_rtaware(mapper);
    else if (mmode == 7) dealloc_mapnode_95_rtaware(mapper);
    else if (mmode == 5) dealloc_mapnode(mapper);
    
    if (rmode == 1) dealloc_dc_lat(router);
    else if (rmode == 4) dealloc_dc_95_lat(router);
    else if (rmode == 2) dealloc_dc_cost(router);
    else dealloc_dc(router);

    // major step 4. Perform computations for comparative benefits table
    
    // case 1. The mmode-rmode of interest.
    
    if (i == 0) {
      printf("Computing statistics for interesting m/r modes\n");
      // pre-compute average costs per client
      for (c=0; c<p.C; c++) {
	
	double vol_c = 0.0;
	for (n=0; n<p.N; n++) vol_c += p.vol[c*N + n];
	
	if (vol_c >= 0.000001) { 
	  // client needs to have sent some traffic
	  // to be included in the analysis
	  avglats[c] = compute_avglats_c(p, alpha, beta, c);
	  avgcost[c] = compute_avgcost_c(p, alpha, beta, c);
	}
      }
    }
    
    // case 2. Other mmode-rmode combinations
    
    else if (i >= 1) {
      
      printf("Comparing statistics for scheme %d\n", i);

      // Zero out statistics arrays so things can just be added later on
      init_counts_diffs(count_cli, count_vol, difflats_cli, diffcost_cli, difflats_vol, diffcost_vol);
      
      // counting non-zero volume clients for averaging later on
      int nonzerovol_cli = 0; 

      // run comparison computation here
      // comparison between (old_alpha, old_beta) and (alpha, beta)
      


      for (c=0; c<p.C; c++) { // for every client
	
	// determine client volume here, because this will be useful everywhere later on
	double vol_c = 0.0;
	for (n=0; n<p.N; n++) vol_c += p.vol[c*N + n];
	
	if (vol_c >= 0.000001) {
	  // client needs to have sent some traffic
	  // to be part of the analysis
	
	  // compute averaged cost and latency metrics
	  avglats_c = compute_avglats_c(p, alpha, beta, c);
	  avgcost_c = compute_avgcost_c(p, alpha, beta, c);
	  

	  // Set an index according to the case where this client falls.
	  int rindex = -1; // initialization
	
	  //printf("client %d avgcost %lf avglats %lf refcost %lf reflats %lf\n", c, avgcost_c, avglats_c, avgcost[c], avglats[c]);
	
	  if (fabs(avgcost_c - avgcost[c]) <= COSTTOL &&
	      fabs(avglats_c - avglats[c]) <= LATSTOL) rindex = 4;

	  else if (avgcost_c > avgcost[c] && avglats_c > avglats[c]) rindex = 0;
	  else if (avgcost_c > avgcost[c] && avglats_c < avglats[c]) rindex = 1;
	  else if (avgcost_c < avgcost[c] && avglats_c > avglats[c]) rindex = 2;	
	  else if (avgcost_c < avgcost[c] && avglats_c < avglats[c]) rindex = 3;
	
	  if (rindex == -1) {
	    printf("Could not classify client %d into any useful class!\n", c);
	    exit(0);
	  }
	
	
	  // stat 1. number of clients in class
	  count_cli[rindex] += 1;
	  nonzerovol_cli += 1;
	
	  // stat 2. clients in this class by volume
	  count_vol[rindex] += vol_c;
	
	  // stat 3. difference in cost by client
	  diffcost_cli[rindex] += fabs(avgcost_c - avgcost[c]);
	
	  // stat 4. difference in cost by volume
	  diffcost_vol[rindex] += (fabs(avgcost_c - avgcost[c]) * vol_c);
	
	  // stat 5. difference in latency by client
	  difflats_cli[rindex] += fabs(avglats_c - avglats[c]);
	
	  // stat 6. difference in latency by volume
	  difflats_vol[rindex] += (fabs(avglats_c - avglats[c]) * vol_c);

	} // analysis for client who sent traffic
	
      } // for every client
      
      // Average the added client metrics here      
      int clas;
      for (clas = 0; clas < 5; clas++) {
	if (count_cli[clas] >= 0.000001) {
	  diffcost_cli[clas] /= count_cli[clas];
	  difflats_cli[clas] /= count_cli[clas];
	}
	// I assume below that p.C and p.TT are necessarily > 0. It's an error otherwise.
	assert(p.TT >= 0.000001 && p.C >= 1);
	count_vol[clas] = 100 * count_vol[clas] / p.TT;
	diffcost_vol[clas] /= p.TT;
	difflats_vol[clas] /= p.TT;
	count_cli[clas] = count_cli[clas] * 100 / nonzerovol_cli;
      }
      
      // Print the averaged difference metrics for this mapping/routing mode
      
      print_diff_stats(mmode, rmode, count_cli, count_vol, diffcost_cli, diffcost_vol, difflats_cli, difflats_vol);
      
    } // comparitive computations for different mapping/routing modes
  } // for each mapping and routing mode compared


 // Finally: deallocate stuff.
  fclose(ftraffic);
  free(alpha);
  free(beta);
  free(int_alpha);
  free(int_beta);
  free(avglats);
  free(avgcost);
  
  MSK_deleteenv(&env);

}

