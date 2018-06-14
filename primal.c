#include <stdio.h>
#include "mosek.h" 
#include "helper.h"
#include "params.h"
#include "mapnode.h"
#include "diag.h"
#include "dc.h"
#include "math.h"
#include "stats.h"

// including feasibility check function here, no header (only one function outside)
int check_feasibility(struct params_t p, double* alpha, double* beta);

// prints MOSEK output to the terminal
static void MSKAPI printstr(void *handle, char str[]) {
  printf("%s", str);
}

// convergence checking function
int check_convergence(MSKrealt objvalue, MSKrealt* old_objvalue, double* beta, double* old_beta, int beta_size) {

  // check objective value differences from one iteration to the next
  // for convergence.

  double change_thresh = 0.01; // absolute threshold for convergence
  int not_converged = 0;

  if (*old_objvalue >= 0.000001) {
    if ((fabs(objvalue - *old_objvalue)/(*old_objvalue)) > 0.01) // 1% threshold for convergence
      not_converged = 1;
  }
  else if (fabs(objvalue - *old_objvalue) >= change_thresh) // absolute threshold for first iter
    not_converged = 1;

  fprintf(stderr, "Old objvalue %lf New objvalue %lf\n", *old_objvalue, objvalue);

  *old_objvalue = objvalue; // for comparison next iteration

  return not_converged;

  // strategy 2. check for convergence based on
  // values of beta between one iteration and the next.
  // But primal optimal variables need not converge.

  /* double beta_change_thresh = 0.0001; */
  
  /* // check for convergence: max difference between old betas and new betas is < a small number */
  /* int i; */
  /* double current_max = 0.0; */
  /* for (i=0; i<beta_size; i++) { */
  /*   if (fabs(beta[i] - old_beta[i]) > current_max) */
  /*     current_max = fabs(beta[i] - old_beta[i]); */
  /*   old_beta[i] = beta[i]; */
  /* } */
  /* if (current_max >= beta_change_thresh) // not converged */
  /*   return 1; */
  /* return 0; // converged case */

}

// takes care of creating a mosek environment, task, and setting up the problem to solve in mosek
int primal(MSKenv_t* env, struct params_t p, struct msk_stats_t* s, int* feasible, double* decomp_objvalue, int* lastiter) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  // initialize optimization variables (primal and dual)
  double* beta = (double*)malloc((I*J*C + I*J)*sizeof(double));
  double* alpha = (double*)malloc((I*C*N + I*J)*sizeof(double));

  // The iterations start with a fixed value of alpha
  // init_beta(p, beta);
  init_alpha(p, alpha); // can depend on which mapping is used. for now, random.
  
  // iteration count for best response game
  int iter;

  // allocate mapnode structures
  struct msk_problem_t* mn = alloc_mapnode(p);
  struct msk_problem_t* dc = alloc_dc(p);

  MSKrescodee r = MSK_RES_OK;
  if (r == MSK_RES_OK) {

    MSKrealt objvalue;
    MSKrealt old_objvalue = 1.0; // initialization to trigger at least two function calls 
                                 // to check_convergence

    // Initialize diagnostics information
    /* int diag_obj_alpha = diag_init("obj_alpha.txt"); */
    /* int diag_obj_beta  = diag_init("obj_beta.txt"); */
    /* int diag_alpha = diag_init("alpha.txt"); */
    /* int diag_beta  = diag_init("beta.txt"); */

    // Write initial beta values, for the record
    // diag_write_list(diag_beta, 0, beta, I*J*C);
    
    // initialize for convergence check function
    // double* old_beta = (double*)calloc(I*J*C, sizeof(double));
    double* old_alpha = (double*)calloc(I*C*N, sizeof(double));

    // set feasibility check variable
    *feasible = 1;
    
    // initialize iteration count
    iter = 0;
    
    //while (check_convergence(objvalue, &old_objvalue, beta, old_beta, I*J*C) > 0 && *feasible && iter < MAXITER) {
    while (check_convergence(objvalue, &old_objvalue, alpha, old_alpha, I*C*N) > 0 && *feasible && iter < MAXITER) {

      *feasible = 0;

      iter ++;

      // **** STUFF WHEN BETA IS INITIALIZED BELOW *** //

      /* //fprintf(stderr, "*****************************************\n"); */
      /* //fprintf(stderr, "*** Primal decomposition iteration %d ***\n", iter); */
      /* //fprintf(stderr, "*****************************************\n"); */

      /* // setup mapnode optimization variables */
      /* setup_mapnode(p, mn, mq, beta); */

      /* // perform optimization using variables in mp structure */
      /* if (optimize_mapnode(&env, p, mn, mq, alpha, &objvalue, /\*iter==(maxiter-1)*\/ 0)) break; */

      /* // write diagnostics for objective value and alpha */
      /* /\* diag_write(diag_obj_alpha, iter, objvalue); *\/ */
      /* /\* diag_write_list(diag_alpha, iter, alpha, I*C*N); *\/ */

      /* // setup dc optimization variables */
      /* setup_dc(p, dc, alpha); */

      /* // perform optimization using variables in md structure */
      /* if (optimize_dc(&env, p, dc, beta, &objvalue, /\*iter==(maxiter-1)*\/ 0)) break; */

      /* // write diagnostics for objective value and beta */
      /* /\* diag_write(diag_obj_beta, iter+1, objvalue); *\/ */
      /* /\* diag_write_list(diag_beta, iter+1, beta, I*J*C); *\/ */

      // *** STUFF WHEN ALPHA IS INITIALIZED *** //

      setup_dc(p, dc, alpha);
      if (optimize_dc(env, p, dc, beta, &objvalue, /*iter==(maxiter-1)*/ 0)) break;
      
      fprintf(stderr, "One round of routing optimization, iteration %d\n", iter);
      
      setup_mapnode(p, mn, beta);
      if (optimize_mapnode(env, p, mn, alpha, &objvalue, /*iter==(maxiter-1)*/ 0)) break;

      fprintf(stderr, "One round of mapping optimization, iteration %d\n", iter);

      // *** COMMON IRRESPECTIVE OF INITIAL SETTING *** //

      // check feasibility test on alphas and betas returned finally
      *feasible = check_feasibility(p, alpha, beta);

      // include data into frequency distributions if feasible
      if (*feasible)
	updatecdf_interim(s, alpha, beta);

    }

    //free(old_beta);
    free(old_alpha);
    
    if (*feasible) {
      updatecdf_finals(s, alpha, beta, iter);
      *decomp_objvalue = objvalue;
    }

    // Set last iteration count
    *lastiter = iter;

    // close diagnostic streams
    /* diag_end(diag_obj_alpha); */
    /* diag_end(diag_obj_beta); */
    /* diag_end(diag_alpha); */
    /* diag_end(diag_beta); */
    
  }

  else {
    /* In case of an error print error code and description. */
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    
    printf("An error occurred while optimizing.\n");
    MSK_getcodedesc(r, symname, desc);
    printf("Error %s - '%s'\n", symname, desc);
  }

  // deallocate mapnode and DC structures
  dealloc_mapnode(mn);
  dealloc_dc(dc);
  
  // free variables
  free(alpha);
  free(beta);
  
  return 0;
  
}
