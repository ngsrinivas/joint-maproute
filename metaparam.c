/*********************************************************************
Generate parameters for a joint mapping+routing optimization problem.
*********************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "params.h"
#include "helper.h"
#include "stats.h"
#include "mosek.h"
#include "global.h"

// prints MOSEK output to the terminal
static void MSKAPI printstr(void *handle, char str[]) {
  printf("%s", str);
}

// Including prototypes of functions called in global.c and primal.c
// Not in a separate header, because these are the only functions from
// these files.
//int global(double* price, double* perf, double* vol, double* w, double* eps, double* cap, int I, int J, int N, int C, double mu, double K, double* global_objvalue);

int primal(MSKenv_t *env, struct params_t p, struct msk_stats_t* s, int* feasible, double* decomp_objvalue, int* lastiter);

/*** Memory allocators ***/

// General double array allocator
double* alloc_double(double** ptr, int size) {
  *ptr = (double*)malloc(size * sizeof(double));
  if (*ptr == 0) { printf("memalloc error!\n"); }
  return *ptr;
}

int alloc_foo(double** price, double** perf, double** vol, double** w, double** eps, double** cap, int I, int J, int N, int C) {
  if(alloc_double(price, I * J) == 0) return -1;
  if(alloc_double(perf, I * J * C) == 0) return -1;
  if(alloc_double(vol, C * N) == 0) return -1;
  if(alloc_double(w, I) == 0) return -1;
  if(alloc_double(eps, I) == 0) return -1;
  if(alloc_double(cap, I * J) == 0) return -1;
  return 0;
}

int dealloc_foo(double* price, double* perf, double* vol, double* w, double* eps, double* cap) {
  free(price);
  free(perf);
  free(vol);
  free(w);
  free(eps);
  free(cap);
  return 0;
}

/**** Error reporting ****/

int err_fun(int code) {
  switch(code) {
  case 1:
    printf("Memory allocation error\n");
    break;
  default:
    printf("Default error\n");
  }
  return 0;
}

/*** Data generators for problem parameters ***/

// General random number generator. Needs bounds for generation
// and number of elements.

int gen_double(double* ptr, int size, double lb, double ub) {
  int i;
  for (i=0; i<size; i++) {
    ptr[i] = lb + ((ub-lb) * rand() / (RAND_MAX + 1.0)); // assume uniform distribution
    //ptr[i] = ub;
  }
}

// generate traffic volumes by assuming one client always resolves from one mapnode
int gen_vols(int C, int N, double* ptr, double lb, double ub) {
  int c, n;
  for (c=0; c<C; c++) {
    // for every client, choose a map node first
    int chosen_mapnode = (N-1) * (double)rand() / (RAND_MAX + 1.0);
    for (n=0; n<N; n++) // start by zeroing all mapnode traffic
      ptr[c*N + n] = 0.0;
    // send random amount of traffic to one chosen mapnode for this client
    ptr[c*N + chosen_mapnode] = lb + ((ub-lb) * rand() / (RAND_MAX + 1.0));
  }
  return 0;
}

int normalize(double* w, int size) {
  double sum = 0.0;
  int i;
  for(i=0; i<size; i++) sum += w[i];
  for(i=0; i<size; i++) w[i] = w[i] / sum;
  return 0;
}

// generate performance values that stress test bad local optima for simulations
int gen_perf(double* perf, int size, double lb, double ub) {
  int i, range;
  double new_lb = ub/1.5;
  double new_ub = lb*1.5;
  for(i=0; i<size; i++) {
    if (((rand() * 1.0)/(RAND_MAX + 1.0)) > 0.5)
      perf[i] = new_lb + ((ub-new_lb) * rand() / (RAND_MAX + 1.0));
    else perf[i] = lb + ((new_ub-lb) * rand() / (RAND_MAX + 1.0));
  }
  return 0;
}

// generate DC relative split in two ranges
int gen_wi(double* w, int I, double lb, double ub) {
  int i;
  double new_lb = ub/1.1;
  double new_ub = lb*1.1;
  for (i=0; i<I; i++) {
    if (((rand() * 1.0)/(RAND_MAX + 1.0)) > 0.5)
      w[i] = new_lb + ((ub-new_lb) * rand() / (RAND_MAX + 1.0));
    else w[i] = lb + ((new_ub-lb) * rand() / (RAND_MAX + 1.0));
  }
  normalize(w, I);
  return 0;
}

int gen_foo(double* price, double* perf, double* vol, double* w, double* eps, double* cap, int I, int J, int N, int C) {
  srand((unsigned int)time(0));

  // choosing price range variation to be fairly small
  // links differ just by performance, to induce bad local optima
  // gen_double(price, I*J, 0.0, 50.0); // price range per 1 request per unit time
  gen_double(price, I*J, 0.5, 1.0);

  // specialized performance values generation to stress test bad local optima
  // gen_double(perf, I*J*C, 50.0, 3000); // in milliseconds
  gen_perf(perf, I*J*C, 50.0, 3000); // in milliseconds. 

  // request rate variation narrowed down
  // gen_double(vol, C*N, 500, 5000); // in request rates
  // gen_vols(C, N, vol, 500, 5000); // in request rates
  gen_vols(C, N, vol, 3000, 5000);

  // DC relative load splits, now into two ranges
  // gen_double(w, I, 1.0, 10.0); // some relative splits, will be normalized
  // normalize(w, I);
  gen_wi(w, I, 1.0, 10.0); // relative splits generation in two ranges now!

  //gen_double(eps, I, 0.0, 0.02); // allowed deviation from target relative split
  gen_double(eps, I, 0.0, 1.0); // allowed deviation from target relative split

  // gen_double(cap, I*J, 1000000, 2000000); // capacities in request rates per second
  gen_double(cap, I*J, 13000000.0, 19000000.0);

  return 0;
}

/****** Output functions ********/

int print_double(FILE* fp, double* ptr, int size) {
  int i;
  fprintf(fp, "%lf", ptr[0]);
  for(i=1; i<size; i++)
    fprintf(fp, " %lf", ptr[i]);
  fprintf(fp, "\n");
}

int output_foo(FILE* fp, double* price, double* perf, double* vol, double* w, double* eps, double* cap, int I, int J, int N, int C, double mu, double K) {
  
  fprintf(fp, "%d %d %d %d\n", I, J, C, N);
  fprintf(fp, "%lf %lf\n", mu, K);
  print_double(fp, price, I*J);
  print_double(fp, perf, I*J*C);
  print_double(fp, vol, C*N);
  print_double(fp, w, I);
  print_double(fp, eps, I);
  print_double(fp, cap, I*J);
  
  return 0;
}

double compute_reldiff(double obj1, double obj2) {
  return (obj2-obj1)*100.0/obj1;
}


// sets link count info for every DC
int set_synth_lnum(int* lnum, int I, int J) {
  int i;
  for (i=0; i<I; i++) lnum[i] = J; // set link count to J for every DC
}
  
// One complete problem iteration: including global problem and primal-dual decomposition
int metaparam_main(MSKenv_t* env, int I, int J, int C, int N, double mu, double K, struct msk_stats_t *s, int problem_index) {
  double *price, *perf, *vol, *w, *eps, *cap;
  double global_objvalue, decomp_objvalue;
  int lastiter, feasible;

  // Allocate memory for parameter vectors
  if(alloc_foo(&price, &perf, &vol, &w, &eps, &cap, I, J, N, C) < 0) return err_fun(1);
  
  // generate parameter vectors one by one -- soon replaced by actual input reading
  gen_foo(price, perf, vol, w, eps, cap, I, J, N, C);

  // write global problem file and get back objective value
  // global(price, perf, vol, w, eps, cap, I, J, N, C, mu, K, &global_objvalue);
  
  // setup structure for primal decomposition gradient
  struct params_t p;
  p.price = price;
  p.perf = perf;
  p.vol = vol;
  p.w = w;
  p.eps = eps;
  p.cap = cap;
  p.I = I;
  p.J = J;
  p.N = N;
  p.C = C;
  p.mu = mu;
  p.K = K;
  p.L = 1.0;
  p.total_links = I*J;
  p.lnum = (int*)malloc(I * sizeof(int));
  set_synth_lnum(p.lnum, p.I, p.J);
  p.TT = total_traffic(vol, C*N);

  fprintf(stderr, "Finished generating and assigning parameters\n");

  // setup and run global optimization
  double* X = (double*)malloc((I*J*C + I*J)*sizeof(double));
  struct msk_problem_t* xp = alloc_global(p);
  setup_global(p, xp);
  if (optimize_global(env, p, xp, X, &global_objvalue, 0)) { // action when optimization fails
    // nothing here for now,
    // but will need to write down the global problem somewhere for 
    // later diagnosis separately.
  }
  // can clear out global optimization memory right here, because the objective
  // is the only thing we need for here.
  free(X);
  dealloc_global(xp);
  fprintf(stderr, "Finished solving global problem\n");
  
  // run primal decomposition iterations -- may remove this under common optimization framework
  primal(env, p, s, &feasible, &decomp_objvalue, &lastiter);

  // Update objective difference distribution if everything is successful
  if (global_objvalue >= 0.000001 && feasible == 1 && lastiter < MAXITER) {
    fprintf(stderr, "Extracted global and decomposed problem results successfully [%lf %lf]\n", global_objvalue, decomp_objvalue);
    printf("SUCCESS\n");
    double objdiff = compute_reldiff(global_objvalue, decomp_objvalue);
    updatecdf_objdiff(s, objdiff);
  }

  // Otherwise, perform error recovery/consolidation actions
  else {

    printf("FAIL\n");

    // case 0. decomposition took too long to converge
    if (lastiter >= MAXITER) {
      fprintf(stderr, "E: Problem %d: Iterations took too long to converge\n", problem_index);
      
      // I'm not quite sure why this would happen, but somtimes I do notice
      // that the problem does not terminate quickly. Putting this in here
      // just in case.
      
      // write problem parameters and global problem into problem specific files
      char command[200];
      sprintf(command, "mv " GLOBAL_PROBLEM_FILENAME " problematic/%d_global.opf", problem_index);
      system(command);

      // write out problem parameters separately
      char filename[100];
      sprintf(filename, "problematic/%d_decomp_divergent.txt", problem_index);
      FILE* fp = fopen(filename, "w");
      if (fp != 0) {
	fprintf(fp, "%d\n", lastiter); // iteration number at which infeasibility was found
	output_foo(fp, price, perf, vol, w, eps, cap, I, J, N, C, mu, K);
	fclose(fp);
      }
      else
	fprintf(stderr, "E: Miserably Failed to open problem parameters file\n");
    }

    // case 1. if the decomposed problem is infeasible
    if (feasible == 0) {
      fprintf(stderr, "E: Problem %d: infeasibility in decomposed problem\n", problem_index);

      // output problem parameters into a problem specific file
      char filename[100];
      sprintf(filename, "problematic/%d_decomp_infeasible.txt", problem_index);
      FILE* fp = fopen(filename, "w");
      if (fp != 0) {
	fprintf(fp, "%d\n", lastiter); // iteration number at which infeasibility was found
	output_foo(fp, price, perf, vol, w, eps, cap, I, J, N, C, mu, K);
	fclose(fp);
      }
      else
	fprintf(stderr, "E: Miserably Failed to open problem parameters file\n");

      // save the two problem description files for inspection too
      char command1[200], command2[200];
      sprintf(command1, "mv problem_alpha.opf problematic/%d_alpha_infeasible.opf", problem_index);
      sprintf(command2, "mv problem_beta.opf  problematic/%d_beta_infeasible.opf", problem_index);
      system(command1);
      system(command2);
    }

    // case 2. extracting global objective was unsuccessful
    if (global_objvalue <= 0.000001) {
      fprintf(stderr, "E: Problem %d: could not extract global objective value successfully!\n", problem_index);

      // copy global problem file into a problem specific file
      char command[200];
      sprintf(command, "mv " GLOBAL_PROBLEM_FILENAME " problematic/%d_global.opf", problem_index);
      system(command);
    }
  }

  // deallocate memory
  dealloc_foo(price, perf, vol, w, eps, cap);
  
  return 0;

}


int main(int argc, char** argv) {

  int I, J, C, N; // DCs, links, clients, mapnodes resp.
  double global_objvalue, decomp_objvalue; // global and decomposed objective values
  int lastiter; // iteration number of the decomposed problem
  int PROBLEM_COUNT; // number of problems to generate
  
  // get number of problems (instances) to run for.
  scanf("%d", &PROBLEM_COUNT);

  // get parameters defining problem size
  scanf("%d", &I);
  scanf("%d", &J);
  scanf("%d", &N);
  scanf("%d", &C);

  // get global constants
  double mu, K;
  scanf("%lf", &mu);
  scanf("%lf", &K);

  // setup structure for CDF statistics computation
  struct msk_stats_t s;
  alloc_cdf(&s, I, J, C, N, PROBLEM_COUNT);
  
  // setup mosek environment
  MSKenv_t env = NULL;
  init_mosek_env(&env);

  printf("Starting %d iterations\n", PROBLEM_COUNT);

  int i, feasible;
  for (i=0; i<PROBLEM_COUNT; i++) {
    sleep(1); // artificial delay between successive problems
    printf("%d ", i);
    fprintf(stderr, "****\n%d\n", i);
    metaparam_main(&env, I, J, C, N, mu, K, &s, i);
    fprintf(stderr, "\n"); // extra line for readability
  }

  // clear mosek environment
  MSK_deleteenv(&env);

  // write CDF statistics; free statistics memory when done
  write_cdfs(&s);
  dealloc_cdf(&s);


  return 0;
}

