#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "params.h"
#include "mosek.h"
#include "inputs.h"
#include "mapnode.h"
#include "dc.h"
#include "result-compute.h"
#include<math.h>
#include<string.h>
#include "cdf.h"

#define MAXITER 100
#define NUMCDFBINS 1000

#define NORMALVAR_TOLERANCE 1e-08
#define DEBUG_RNG 0
#define DEBUG_CDF 0

unsigned int seed1, seed2;

// Checking convergence using objective values.
int check_convergence(double new_obj, double* old_obj) {

  int not_converged = 1;
  
  if (*old_obj >= 0.000001) {
    if ((fabs(new_obj - *old_obj)/(*old_obj)) < 0.01) not_converged = 0;
  }

  *old_obj = new_obj;
  
  return not_converged;
}

// Record original traffic volumes, these are the means
// for independent price generation later on.
int record_original_prices(struct params_t p, double* prices) {
  
  int I = p.I;
  int J = p.J;
  
  memcpy((void*)prices, (void*)p.price, I*J*sizeof(double));
  
  return 0;
}

int get_01normal(double* z0, double* z1) {

  double u, v, s;
  
  do {
    u = -1.0 + (rand_r(&seed1) * 2.0 / RAND_MAX);
    v = -1.0 + (rand_r(&seed2) * 2.0 / RAND_MAX);
    
    //if (DEBUG_RNG) printf("u: %lf v:%lf\n", u, v);

    s = u*u + v*v;

  } while (s <= 0.000001 || s >= 1.000000);
  
  double f = sqrt(-2 * log(s) / s);
  *z0 = u * f;
  *z1 = v * f;
  
  return 0;

}

// returns a clipped gaussian random number (clipped on
// the positive half quadrant) with the provided mean
// and variance
double next_gaussian(double mean, double variance) {
  
  double z0, z1;
  double normalvar = -1.0;
  
  int itercount = 0;

  do {
    get_01normal(&z0, &z1);
    normalvar = mean + variance * z0;
    if (normalvar <= NORMALVAR_TOLERANCE) normalvar = mean + variance * z1;
    
    if (DEBUG_RNG) printf("normalvar: %e z0: %lf z1: %lf condition: %d\n", normalvar, z0, z1, normalvar <= NORMALVAR_TOLERANCE);
    
    itercount ++;
  } while (normalvar <= NORMALVAR_TOLERANCE && itercount < 100);
  
  if (itercount >= 100) normalvar = mean;
  
  printf("normalvar %lf mean %lf\n", normalvar, mean);
  
  return normalvar;
  
}


// initialize seeds for independent RNGs
int init_seeds() {
  
  seed1 = (unsigned)time(0);

  FILE* fp = fopen("/dev/random", "r");
  //fscanf(fp, "%d", &seed1);
  fscanf(fp, "%d", &seed2);
  fclose(fp);
  
  return 0;

}

// Use the variance and generate prices *independently* for the
// links. Use the mean from the `prices` argument. The variance is
// a fraction of the mean for the particular link, and the resulting
// random variable is capped on the bottom by 0.0 of course. Samples
// are independent.
int set_prices(struct params_t* p, double* prices, double varfrac) {
  
  int I = p->I;
  int J = p->J;
  int i, j;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p->cap[i*J + j] >= 0.000001) { // valid link
	if (DEBUG_RNG) printf("mean: %lf variance: %lf \n", prices[i*J + j], varfrac * prices[i*J + j]);
	p->price[i*J + j] = next_gaussian(prices[i*J + j], varfrac * prices[i*J + j]);
      }
    }
  }
  
  return 0;

}

// Do map-route iterations or one-time optimizations based on input modes.
int main(int argc, char** argv) {
  
  // attempt to get to a pareto optimal curve.
  struct params_t p;
  struct cdf_t *wperf, *objperf, *bprice;
  
  if (argc < 15) {
    printf("Usage: %s basic-inputs linknum-info traffic-file perf-file caps-file pricing-file mapping-mode[1-5] routing-mode[1-3] outfile_wperf outfile_objperf outfile_tprice outfile_gen varfrac nsamples\n", argv[0]);
    exit(0);
  }

  init_seeds();
  
  // process inputs and keep a few things handy
  input_master_single_instance(&argv[1], &p);

  int mmode = atoi(argv[7]); // assumed (and NOT CHECKED) to be [1-6]
  int rmode = atoi(argv[8]); // assumed (and NOT CHECKED) to be [1-3]
  
  // set up file descriptors for CDF
  FILE* fwperf = fopen(argv[9], "w");
  FILE* fobjperf = fopen(argv[10], "w");
  FILE* fbprice = fopen(argv[11], "w");
  FILE* fgen = fopen(argv[12], "a");
  if (fwperf == NULL) {
    printf("Could not open output file %s !\n", argv[9]);
    exit(0);
  }
  if (fobjperf == NULL) {
    printf("Could not open output file %s !\n", argv[10]);
    exit(0);
  }
  if (fbprice == NULL) {
    printf("Could not open output file %s !\n", argv[11]);
    exit(0);
  }			
  if (fgen == NULL) {
    printf("Could not open output file %s !\n", argv[12]);
    exit(0);
  }
  
  // Setup sampling information
  double varfrac = atof(argv[13]);
  int nsamples = atoi(argv[14]); // number of samples overall
  
  
  // Now, turning to business...
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
  
  // allocate cdf structures
  wperf = cdf_init_doublecdf(nsamples, NUMCDFBINS);
  objperf = cdf_init_doublecdf(nsamples, NUMCDFBINS);
  bprice = cdf_init_doublecdf(nsamples, NUMCDFBINS);
  
  // mosek problem setup related structures
  struct msk_problem_t *mapper, *router;
  double objvalue, old_objvalue=0.0;
  double mapping_objvalue;

  // Record original prices somewhere; we'll be
  // using this to determine the prices later on.
  double* prices = (double*)malloc(I*J*sizeof(double));
  record_original_prices(p, prices);
  
  if (rmode == 4) {
    fprintf(stderr, "Need to fix initializations first\n");
    exit(0);
  }
  
  if (mmode == 4 || mmode == 5 || mmode == 7) {
    
    init_alpha(p, alpha); // round robin initialization

  }
  
  // Mapping memory allocations

  if(mmode == 3) mapper = alloc_mapnode_capped_closest(p);
  else if(mmode == 6) mapper = alloc_mapnode_95_capped_closest(p);
  else if (mmode == 4) {
    mapper = alloc_mapnode_rtaware(p);
  }  
  else if(mmode == 5) {
    mapper = alloc_mapnode(p);
  }
  else if (mmode == 7) {
    mapper = alloc_mapnode_95_rtaware(p);
  }


  // Routing memory allocations
  
  if (rmode == 1) router = alloc_dc_lat(p);
  else if (rmode == 4) router = alloc_dc_95_lat(p);
  else if (rmode == 2) router = alloc_dc_cost(p);
  else router = alloc_dc(p);

  // loop through samples

  int sample;

  for (sample = 0; sample < nsamples; sample++) { // for each sample
  
    printf("Sample=%d \n", sample);
    
    // (major) step 1. Set up prices based on the provided variance factor
    set_prices(&p, prices, varfrac);

    // Any optimization steps are needed only if the routing mode is
    // not 1, or if the sample index is 0. So, if the routing mode is
    // performance only, optimal mapping and routing variables are
    // computed only for the first sample. Only prices are recomputed
    // for the later samples. Saves a TON of time.
    
    if ((rmode != 1 && rmode != 4) || sample == 0) {

      // (major) step 2. Do one time mapping.

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

  
      // (major) step 3. Perform one-time routing optimization
      // or best response depending on the mapping mode used.
    
      // one time routing

      if (mmode <= 3 || mmode == 6) {
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

      // Best response iterative routing
    
      else {
	int numiter = 0;
	old_objvalue = 0.0;

	// debugging code to skip optimization steps
	if (DEBUG_CDF) {
	  old_objvalue = rand_r(&seed1) / (RAND_MAX + 1.0);
	  objvalue = old_objvalue;
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

    } // end optimization-required block (rmode != 1 || sample == 0)
    
    // (major) step 4. record stuff into generic output file
    
    double objprice = compute_total_price_alphabeta(&p, alpha, beta);
    double objperfval = 0.0;
    
    if (rmode == 1) objperfval = objvalue / p.K;
    else if (rmode == 4 || rmode == 2) objperfval = compute_objperf_alphabeta(&p, alpha, beta);
    else if (rmode == 3) objperfval = (objvalue - objprice) / p.K;
    else {
      printf("Can't infer objperf from routing objective value; implement this function for result computation!\n");
      exit(0);
    }
    
    double wperfval = compute_perf_alphabeta(&p, alpha, beta);
    double TT = p.TT;
    double bpriceval = objprice * 1000000 / (65.8 * TT);
    double wrtt = compute_just_wrtt_alphabeta(&p, alpha, beta);
    
    if (DEBUG_CDF) {
      bpriceval = rand_r(&seed1);
      objperfval = rand_r(&seed1);
      wperfval = rand_r(&seed1);
    }

    // data format (7 fields)
    // sample_index objperf(ms) wperf(ms) wrtt(ms) tprice($/hour) bprice ($/GB) routing_objvalue_value
    fprintf(fgen, "%d %lf %lf %lf %lf %lf %lf\n", sample, objperfval, wperfval, wrtt, 3600*objprice, bpriceval, objvalue);
    
    // (major) step 5. add items into CDFs
    cdf_additem(objperf, objperfval);
    cdf_additem(wperf, wperfval);
    cdf_additem(bprice, bpriceval);

  } // for each sample
    
  // Finally: deallocate stuff.
  free(alpha);
  free(beta);

  // write cdf and free cdf memory
  cdf_compute_doublecdf(wperf);
  cdf_compute_doublecdf(objperf);
  cdf_compute_doublecdf(bprice);
  cdf_writecdf(fwperf, wperf);
  cdf_writecdf(fbprice, bprice);
  cdf_writecdf(fobjperf, objperf);
  cdf_dealloc_doublecdf(wperf);
  cdf_dealloc_doublecdf(bprice);
  cdf_dealloc_doublecdf(objperf);

  // close file descriptors
  fclose(fwperf);
  fclose(fbprice);
  fclose(fobjperf);
  fclose(fgen);

  if (mmode == 3) dealloc_mapnode_capped_closest(mapper);
  else if (mmode == 6) dealloc_mapnode_95_capped_closest(mapper);
  else if (mmode == 4) dealloc_mapnode_rtaware(mapper);
  else if (mmode == 7) dealloc_mapnode_95_rtaware(mapper);
  else if (mmode == 5) dealloc_mapnode(mapper);
  
  if (rmode == 1) dealloc_dc_lat(router);
  else if (rmode == 4) dealloc_dc_95_lat(router);
  else if (rmode == 2) dealloc_dc_cost(router);
  else dealloc_dc(router);
  
  MSK_deleteenv(&env);

}

