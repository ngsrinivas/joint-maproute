/*** Functions to prepare objective, coefficients, etc. for the global problem ***/

#include<stdlib.h>
#include "params.h"
#include "mosek.h"
#include "helper.h"

#define MY_EPSILON 0.000001

// extern definitions for link penalty term-related coefficients
extern int y_coeff[6];
extern double cap_coeff[6];

// prints MOSEK output to the terminal
static void MSKAPI printstr(void *handle, char str[]) {
  printf("%s", str);
}

// coefficients in the objective for global problem
int prepare_objective_gl(struct params_t p, double* obj_coeff) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double K = p.K;
  double L = p.L;
  double mu = p.mu;
  double TT = p.TT;

  int i, j, c, n;
  double vol_c;

  // X_ijc variables
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	vol_c = 0.0;
	for (n=0; n<N; n++) // because sum_n vol_cn = vol_c, which we need.
	  vol_c += p.vol[c*N + n];
	obj_coeff[i*J*C + j*C + c] = ((L * vol_c * p.price[i*J + j]) + (K * (vol_c/TT) * p.perf[i*J*C + j*C + c]));
      }
    }
  }
  
  // phi_ij variables
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) { 
      if (p.cap[i*J + j] >= 0.000001) {
	obj_coeff[I*J*C + i*J + j] = K * SC_PHI / (p.cap[i*J + j] * p.total_links);
      }
      else {
	obj_coeff[I*J*C + i*J + j] = 0.0;
      }
    }
  }
  
  return 0;
}

// constraint coefficients, set row by row
int prepare_constraints_gl(struct params_t p, int* aptrb, int* aptre, int* asub, double* aval) {

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
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) { // for each coefficient here.
	asub[subindex] = i*J*C + j*C + c;
	aval[subindex] = 0.0;
	for (n=0; n<N; n++)
	  aval[subindex] += p.vol[c*N + n]/TT;
	subindex ++;
      }
    } // for each subscript
    aptre[row] = subindex;
    row ++;
  }

  // set 2. link load penalty based constraints
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) { // for each link
      
      // 6 constraints per link
      int k;
      for (k=0; k<6; k++) { // for each constraint
	aptrb[row] = subindex;
	for (c=0; c<C; c++) { // for each coefficient of X_ijc
	  asub[subindex] = i*J*C + j*C + c;
	  aval[subindex] = 0.0;
	  for (n=0; n<N; n++)
	    aval[subindex] += p.vol[c*N + n];
	  aval[subindex] *= y_coeff[k];

	  subindex ++;
	} // for each coefficient of X_ijc
	
	// Take care of phi_ij coefficient
	asub[subindex] = I*J*C + i*J + j;
	aval[subindex] = -1.0;
	subindex ++;

	// finished k'th constraint of ij'th link
	aptre[row] = subindex;
	row ++;
	
	
      } // for each link constraint (innermost level)
    } // for each link constraint
  } // for each link constraint (outermost level)

  // set 3. Traffic inclusion
  for (c=0; c<C; c++) { // for each constraint
    aptrb[row] = subindex;
    for (i=0; i<I; i++) {
      for (j=0; j<J; j++) { // for each X_ijc coefficient
	asub[subindex] = i*J*C + j*C + c;
	aval[subindex] = 1.0;
	subindex++;
      }
    } // for each coefficient
    aptre[row] = subindex;
    row ++;
  } // for each constraint
	
  //printf("Global: Number of constraints %d %d Number of nonzero coeffs %d %d\n", row, I+6*I*J+C, subindex, 8*I*J*C + 6*I*J);

  return 0;

}

// constraint bounds, set row by row
int prepare_constbounds_gl(struct params_t p, MSKboundkeye* bkc, double* blc, double* buc) {
  
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  int row = 0;

  int i, j, c, n;
  
  // set 1. Load balancing constraints
  for (i=0; i<I; i++) {
    bkc[row] = MSK_BK_RA;
    blc[row] = p.w[i] - p.eps[i];
    buc[row] = p.w[i] + p.eps[i];
    row ++;
  }

  // set 2. Link capacity constraints

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      
      // set up 6 constraint bounds
      int k;
      for (k=0; k<6; k++) {
	bkc[row] = MSK_BK_UP;
	blc[row] = -MSK_INFINITY;
	buc[row] = cap_coeff[k] * p.cap[i*J + j];
	row ++;
      } // k'th constraint for ij'th link
    }
  } // for ij'th link

  // set 3. Traffic inclusion
  for (c=0; c<C; c++) {
    bkc[row] = MSK_BK_FX;
    blc[row] = 1;
    buc[row] = 1;
    row ++; 
  }

  return 0;

}

// Variable bounds
int prepare_varbounds_gl(struct params_t p, MSKboundkeye* bkx, double* blx, double* bux) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  int* lnum = p.lnum;

  int i, j, c, n;

  int var = 0;

  // X_ijc variables
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	if (j < lnum[i]) {
	  if (p.perf[i*J*C + j*C + c] >= 0.000001) {
	    bkx[var] = MSK_BK_LO;
	    blx[var] = 0.0;
	    bux[var] = MSK_INFINITY;
	  }
	  else {
	    bkx[var] = MSK_BK_FX;
	    blx[var] = 0.0;
	    bux[var] = 0.0;
	  }
	}
	else {
	  bkx[var] = MSK_BK_FX;
	  blx[var] = 0.0;
	  bux[var] = 0.0;
	}
	var ++;
      }
    }
  }

  // phi_ij variables
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      bkx[var] = MSK_BK_LO;
      blx[var] = 0.0;
      bux[var] = MSK_INFINITY;
      var ++;
    }
  }

  return 0;

}


// Setup parameters for the global problem
int setup_global(struct params_t p, struct msk_problem_t *mp) {

  // prepare problem related arrays:
  // 1. objective coefficients
  // 2. constraint matrix sparse description arrays
  // 3. constraint bound arrays
  // 4. update provided structure with problem information

  prepare_objective_gl(p, mp->obj_coeff);
  prepare_constraints_gl(p, mp->aptrb, mp->aptre, mp->asub, mp->aval);
  prepare_constbounds_gl(p, mp->bkc, mp->blc, mp->buc);
  prepare_varbounds_gl(p, mp->bkx, mp->blx, mp->bux);

  return 0;

}

// get memory allocated for various structures of coefficients
struct msk_problem_t* alloc_global(struct params_t p) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  // coefficients in objective functions...
  // corresponds to number of variables involved
  double* obj_coeff = (double*)malloc((I*J*C + I*J)*sizeof(double));

  // constraint
  int* aptrb = (int*)malloc((I  + 6*I*J + C) * sizeof(int)); // number of rows
  int* aptre = (int*)malloc((I  + 6*I*J + C) * sizeof(int)); // number of rows
  int* asub  = (int*)malloc((8*I*J*C + 6*I*J) * sizeof(int)); // number of nonzero terms
  double* aval = (double*)malloc((8*I*J*C + 6*I*J) * sizeof(double)); // number of nonzero terms

  // cnostraint bounds
  MSKboundkeye* bkc = (MSKboundkeye*)malloc((I+6*I*J+C)*sizeof(MSKboundkeye));
  double* blc = (double*)malloc((I+6*I*J+C)*sizeof(double));
  double* buc = (double*)malloc((I+6*I*J+C)*sizeof(double));

  // variable bounds
  MSKboundkeye* bkx = (MSKboundkeye*)malloc((I*J*C + I*J)*sizeof(MSKboundkeye));
  double* blx = (double*)malloc((I*J*C + I*J)*sizeof(double));
  double* bux = (double*)malloc((I*J*C + I*J)*sizeof(double));

  struct msk_problem_t *mp = (struct msk_problem_t*)malloc(sizeof(struct msk_problem_t));

  // Memory alloc checking
  if(obj_coeff == 0 || aptrb == 0 || aptre == 0 || asub == 0 || aval == 0 || bkc == 0 || blc == 0 || buc == 0 || 
       bkx == 0 || blx == 0 || bux == 0 || mp == 0) {
    printf("Memory allocation error msk_problem_t in global problem!\n");
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
int dealloc_global(struct msk_problem_t *mp) {
  
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

int optimize_global(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* X, MSKrealt* objvalue, int print_sol) {
  
  MSKtask_t    task = NULL;
  MSKidxt i,j;
  MSKrescodee  r;
  
  int numvar = p.I * p.J * p.C + p.I * p.J;
  int numcon = p.I + 6 * p.I * p.J + p.C;
  int numanz = 8 * p.I * p.J * p.C + 6 * p.I * p.J;
  
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
			     X);
	
	
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
	printf("Global problem: Primal or dual infeasibility certificate found.\n");
	// write problem data to a file for future inspection
	MSK_writedata(task, "problem_X.opf");
	MSK_deletetask(&task);
	rvalue = 1;
	break;
        
      case MSK_SOL_STA_UNKNOWN:
	printf("Global problem: The status of the solution could not be determined: writing current solutions in\n");
	// write problem data to a file for future inspection, if needed
	MSK_writedata(task, "problem_X.opf");
	
	// get complete primal solution
	MSK_getsolutionslice(task,
			     MSK_SOL_ITR,    /* Request the basic solution. */
			     MSK_SOL_ITEM_XX,/* Which part of solution.     */
			     0,              /* Index of first variable.    */
			     numvar,         /* Index of last variable+1.   */
			     X);

	rvalue = 1;
	break;

      default:
	printf("Global problem: Other solution status.");
	// write problem data to a file for future inspection, if needed
	MSK_writedata(task, "problem_X.opf");
	rvalue = 1;
	break;
      }
    }
    else {
      printf("Global problem: Error while optimizing.\n");
      // write problem data to a file for future inspection, if needed
      MSK_writedata(task, "problem_X.opf");
      rvalue = 1;
    }
  }
  
  if (r != MSK_RES_OK) {
    /* In case of an error print error code and description. */      
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    
    printf("Global problem: An error occurred while optimizing.\n");
    MSK_getcodedesc (r,
		     symname,
		     desc);
    printf("Error %s - '%s'\n",symname,desc);

    rvalue = 1;

  }
  
  MSK_deletetask(&task);

  return rvalue;

}

