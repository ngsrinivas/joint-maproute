#include<stdlib.h>
#include "params.h"
#include "mosek.h"
#include "helper.h"

#define MY_EPSILON 0.000001

// prints MOSEK output to the terminal
static void MSKAPI printstr(void *handle, char str[]) {
  printf("%s", str);
}

// coefficients in the objective for mapping node subproblem
int prepare_objective_mn_capped_closest(struct params_t p, double* obj_coeff) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double K = p.K;
  double L = p.L;
  double mu = p.mu;
  double TT = p.TT;

  int i, j, c, n;
  double coeff_icn, curr_perf;

  for (i=0; i<I; i++) {
    for(c=0; c<C; c++) {
      
      // determine best case latency for DC i
      double minlat = MSK_INFINITY;
      for (j=0; j<J; j++) {
	curr_perf = p.perf[i*J*C + j*C + c];
	if (curr_perf >= 0.000001 && curr_perf < minlat) minlat = curr_perf;
      }
	
      // set alpha_icn coefficients here
      for(n=0; n<N; n++)
	obj_coeff[i*C*N + c*N + n] = (p.vol[c*N + n]/TT) * minlat;
    }
  }

  return 0;
}

// constraint coefficients, set row by row
int prepare_constraints_mn_capped_closest(struct params_t p, int* aptrb, int* aptre, int* asub, double* aval) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double TT = p.TT;
  
  int row = 0; // row count in the constraint matrix
  int subindex = 0; // subscript index (same as value index, in the sparse representation)

  int i, j, c, n;
  
  // set 1. load balancing constraints
  for (i=0; i<I; i++) { // for each constraint
    aptrb[row] = subindex;
    for(c=0; c<C; c++) {
      for(n=0; n<N; n++) {
	asub[subindex] = i*C*N + c*N + n;
	aval[subindex] = p.vol[c*N+n];
	subindex ++;
      }
    }
    aptre[row] = subindex;
    row ++;
  }

  // set 2. Traffic inclusion from mapping nodes
  for (c=0; c<C; c++) {
    for (n=0; n<N; n++) { // for each constraint
      aptrb[row] = subindex;
      for (i=0; i<I; i++) {
	asub[subindex] = i*C*N + c*N + n;
	aval[subindex] = 1;
	subindex ++;
      }
      aptre[row] = subindex;
      row ++;
    }
  }

  //printf("Mapping: Number of constraints %d %d Number of nonzero coeffs %d %d\n", row, I+6*I*J+C*N, subindex, 2*I*C*N + 6*I*J*C*N + 6*I*J);

  return 0;

}

// constraint bounds, set row by row
int prepare_constbounds_mn_capped_closest(struct params_t p, MSKboundkeye* bkc, double* blc, double* buc) {
  
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  int row = 0;

  int i, j, c, n;
  
  // set 1. Load balancing constraints
  for (i=0; i<I; i++) {
    // set to sum of link capacities for the link
    double tcap = 0.0;
    for (j=0; j<J; j++) tcap += p.cap[i*J + j];
    bkc[row] = MSK_BK_UP;
    blc[row] = 0.0;
    buc[row] = tcap;
    row ++;
  }

  // set 2. Traffic inclusion
  for (c=0; c<C; c++) {
    for (n=0; n<N; n++) {
      bkc[row] = MSK_BK_FX;
      blc[row] = 1;
      buc[row] = 1;
      row ++;
    }
  }

  return 0;

}

// Variable bounds
int prepare_varbounds_mn_capped_closest(struct params_t p, MSKboundkeye* bkx, double* blx, double* bux) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  int i, j, c, n;

  int var = 0;

  // alpha variables
  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) {
      for (n=0; n<N; n++) {
	bkx[var] = MSK_BK_LO;
	blx[var] = MY_EPSILON; // a small amount of traffic is _always_ sent
	bux[var] = MSK_INFINITY;
	var ++;
      }
    }
  }

  return 0;

}


// Setup parameters for the mapping node subproblem
int setup_mapnode_capped_closest(struct params_t p, struct msk_problem_t *mp) {

  // prepare problem related arrays:
  // 1. objective coefficients
  // 2. constraint matrix sparse description arrays
  // 3. constraint bound arrays
  // 4. update provided structure with problem information

  prepare_objective_mn_capped_closest(p, mp->obj_coeff);
  prepare_constraints_mn_capped_closest(p, mp->aptrb, mp->aptre, mp->asub, mp->aval);
  prepare_constbounds_mn_capped_closest(p, mp->bkc, mp->blc, mp->buc);
  prepare_varbounds_mn_capped_closest(p, mp->bkx, mp->blx, mp->bux);

  return 0;

}

// get memory allocated for various structures of coefficients
struct msk_problem_t* alloc_mapnode_capped_closest(struct params_t p) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  // coefficients in objective functions...
  // corresponds to number of variables involved
  double* obj_coeff = (double*)malloc((I*C*N)*sizeof(double));

  // constraint
  int* aptrb = (int*)malloc((I + C*N) * sizeof(int)); // number of rows
  int* aptre = (int*)malloc((I + C*N) * sizeof(int)); // number of rows
  int* asub  = (int*)malloc((2*I*C*N) * sizeof(int)); // number of nonzero terms
  double* aval = (double*)malloc((2*I*C*N) * sizeof(double)); // number of nonzero terms

  // cnostraint bounds
  MSKboundkeye* bkc = (MSKboundkeye*)malloc((I+C*N)*sizeof(MSKboundkeye));
  double* blc = (double*)malloc((I+C*N)*sizeof(double));
  double* buc = (double*)malloc((I+C*N)*sizeof(double));

  // variable bounds
  MSKboundkeye* bkx = (MSKboundkeye*)malloc((I*C*N)*sizeof(MSKboundkeye));
  double* blx = (double*)malloc((I*C*N)*sizeof(double));
  double* bux = (double*)malloc((I*C*N)*sizeof(double));

  struct msk_problem_t *mp = (struct msk_problem_t*)malloc(sizeof(struct msk_problem_t));

  // Memory alloc checking
  if(obj_coeff == 0 || aptrb == 0 || aptre == 0 || asub == 0 || aval == 0 || bkc == 0 || blc == 0 || buc == 0 || 
       bkx == 0 || blx == 0 || bux == 0 || mp == 0) {
    printf("Memory allocation error msk_problem_t in mapping node subproblem!\n");
    exit(0);
  }


  mp->obj_coeff = obj_coeff;
  mp->aptrb = aptrb;
  mp->aptre = aptre;
  mp->asub = asub;
  mp->aval = aval;
  mp->bkc = bkc;
  mp->blc = blc;
  mp->buc = buc;
  mp->bkx = bkx;
  mp->blx = blx;
  mp->bux = bux;

  return mp;

}

// deallocate memory within problem structure
int dealloc_mapnode_capped_closest(struct msk_problem_t *mp) {
  
  free(mp->obj_coeff);
  free(mp->aptrb);
  free(mp->aptre);
  free(mp->asub);
  free(mp->aval);
  free(mp->bkc);
  free(mp->blc);
  free(mp->buc);
  free(mp->bkx);
  free(mp->blx);
  free(mp->bux);
  free(mp);

  return 0;

}

// Makes mosek api calls for solving the optimization problem, and
// returns results in alpha and objvalue (primal objective 
// value). Mosek is allowed to print stuff only if print_sol is 1.

int optimize_mapnode_capped_closest(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* alpha, MSKrealt* objvalue, int print_sol) {

  MSKtask_t    task = NULL;
  MSKidxt i,j;
  MSKrescodee  r;
  
  int numvar = p.I * p.C * p.N;
  int numcon = p.I + p.C * p.N;
  int numanz = 2 * p.I * p.C * p.N;
  
  int rvalue = 0; // default return value
  
  // Create the optimization task.
  r = MSK_maketask(*env, numcon, numvar, &task);

  // Directs the log task stream to the 'printstr' function.
  if (print_sol > 0)
    if (r == MSK_RES_OK)
      MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
  
  // Give MOSEK an estimate of the size of the input data.
  if (r == MSK_RES_OK)
    r = MSK_putmaxnumvar(task, numvar);
  if (r == MSK_RES_OK)
    r = MSK_putmaxnumcon(task, numcon);
  if (r == MSK_RES_OK)
    r = MSK_putmaxnumanz(task, numanz);

  // Append empty variables and constraints
  if (r == MSK_RES_OK)
    r = MSK_append(task, MSK_ACC_VAR, numvar);
  if (r == MSK_RES_OK)
    r = MSK_append(task, MSK_ACC_CON, numcon);

  for(j=0; j<numvar && r == MSK_RES_OK; ++j) {
    // Set the linear term c_j in the objective.
    if(r == MSK_RES_OK)
      r = MSK_putcj(task, j, mp->obj_coeff[j]);
    
    // Set the bounds on variable j.
    // blx[j] <= x_j <= bux[j]
    if(r == MSK_RES_OK)
      r = MSK_putbound(task,
		       MSK_ACC_VAR, /* Put bounds on variables.*/
		       j,           /* Index of variable.*/
		       mp->bkx[j],      /* Bound key.*/
		       mp->blx[j],      /* Numerical value of lower bound.*/
		       mp->bux[j]);     /* Numerical value of upper bound.*/
  } // for every variable
  
  // Set the bounds on constraints.
  //   for i=1, ..., NUMCON : blc[i] <= constraint i <= buc[i]

  for(i=0; i<numcon && r==MSK_RES_OK; ++i) {  
    // Set constraint bound
    r = MSK_putbound(task,
		     MSK_ACC_CON, /* Put bounds on constraints.*/
		     i,           /* Index of constraint.*/
		     mp->bkc[i],   /* Bound key.*/
		     mp->blc[i],   /* Numerical value of lower bound.*/
		     mp->buc[i]);  /* Numerical value of upper bound.*/
    
    /* Input column j of A */   
    if(r == MSK_RES_OK)
      r = MSK_putavec(task,
		      MSK_ACC_CON,             /* Input row of A.*/
		      i,                       /* Row index.*/
		      mp->aptre[i]-mp->aptrb[i], /* Number of non-zeros in row i.*/
		      mp->asub+mp->aptrb[i],     /* Pointer to column indexes of row i.*/
		      mp->aval+mp->aptrb[i]);    /* Pointer to Values of row i.*/      
  } // for every constraint

  /* Maximize objective function. */
  if (r == MSK_RES_OK)
    r = MSK_putobjsense(task,
			MSK_OBJECTIVE_SENSE_MINIMIZE);

  // All data has been input in the task. 
  // Solve the problem now.
  if (r == MSK_RES_OK) {
    MSKrescodee trmcode;
    
    /* Run optimizer */
    r = MSK_optimizetrm(task, &trmcode);

    /* Print a summary containing information
       about the solution for debugging purposes. */
    MSK_solutionsummary(task, MSK_STREAM_LOG);
    
    if (r == MSK_RES_OK) {
      MSKsolstae solsta;
      int j;
      MSK_getsolutionstatus (task,
			     MSK_SOL_ITR,
			     NULL,
			     &solsta);
      switch(solsta) {
      case MSK_SOL_STA_OPTIMAL:   
      case MSK_SOL_STA_NEAR_OPTIMAL:

	// get complete primal solution
	MSK_getsolutionslice(task,
			     MSK_SOL_ITR,    /* Request the basic solution. */
			     MSK_SOL_ITEM_XX,/* Which part of solution.     */
			     0,              /* Index of first variable.    */
			     numvar,         /* Index of last variable+1.   */
			     alpha);
	
	
	if (print_sol > 0) {
	  printf("Optimal primal solution\n");
	  //for(j=0; j<numvar; ++j)
	  //  printf("x[%d]: %e\n", j, alpha[j]);
	}

	MSK_getprimalobj(task, MSK_SOL_ITR, objvalue);
	
	break;

      case MSK_SOL_STA_DUAL_INFEAS_CER:
      case MSK_SOL_STA_PRIM_INFEAS_CER:
      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
	printf("Mapping problem: Primal or dual infeasibility certificate found.\n");
	// write problem data to a file for future inspection
	MSK_writedata(task, "problem_alpha.opf");
	MSK_deletetask(&task);
	rvalue = 1;
	break;
        
      case MSK_SOL_STA_UNKNOWN:
	printf("Mapping problem: The status of the solution could not be determined: writing current solutions in\n");
	// write problem data to a file for future inspection, if needed
	MSK_writedata(task, "problem_alpha.opf");
	
	// get complete primal solution
	MSK_getsolutionslice(task,
			     MSK_SOL_ITR,    /* Request the basic solution. */
			     MSK_SOL_ITEM_XX,/* Which part of solution.     */
			     0,              /* Index of first variable.    */
			     numvar,         /* Index of last variable+1.   */
			     alpha);

	rvalue = 1;
	break;

      default:
	printf("Mapping problem: Other solution status.");
	// write problem data to a file for future inspection, if needed
	MSK_writedata(task, "problem_alpha.opf");
	rvalue = 1;
	break;
      }
    }
    else {
      printf("Mapping problem: Error while optimizing.\n");
      // write problem data to a file for future inspection, if needed
      MSK_writedata(task, "problem_alpha.opf");
      rvalue = 1;
    }
  }
  
  if (r != MSK_RES_OK) {
    /* In case of an error print error code and description. */      
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    
    printf("Mapping problem: An error occurred while optimizing.\n");
    MSK_getcodedesc (r,
		     symname,
		     desc);
    printf("Error %s - '%s'\n",symname,desc);

    rvalue = 1;

  }
  
  MSK_deletetask(&task);

  return rvalue;

}

