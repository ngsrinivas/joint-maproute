#include "mosek.h"
#include "params.h"
#include "helper.h"

// prints MOSEK output to the terminal
static void MSKAPI printstr(void *handle, char str[]) {
  printf("%s", str);
}

// prepare objective function
int prepare_objective_dc(struct params_t p, double* alpha, double* obj_coeff) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double K = p.K;
  double mu = p.mu;
  double TT = p.TT;
  double L = p.L;

  int i, j, c, n;
  double coeff_ijc;

  // coefficients on beta's
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	coeff_ijc = 0.0;
	for (n=0; n<N; n++) {
	  coeff_ijc += ((L * alpha[i*C*N+c*N+n] * p.price[i*J+j] * p.vol[c*N+n]) + (K * alpha[i*C*N+c*N+n] * (p.vol[c*N+n]/TT) * p.perf[i*J*C+j*C+c]));
	}
	obj_coeff[i*J*C+j*C+c] = coeff_ijc;
      }
    }
  }

  // coefficients for link load penalty variables, one per link. We use 1/total traffic to 
  // average out over the entire traffic, just like we do for weighted RTT.
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

// prepare constraints
int prepare_constraints_dc(struct params_t p, double* alpha, int* aptrb, int* aptre, int* asub, double* aval) {
  
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double TT = p.TT;

  int row = 0; // row count in the constraint matrix
  int subindex = 0; // subscript index (same as value index, in the sparse representation)

  int i, j, c, n;

  // 6 link capacity constraints for every link, in terms of penalty variable.
  // Use Joe's function here directly
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) { // for each link
      
      // In constraints below, we denote link load as y_ij and the 
      // corresponding link utilization penalty variable as phi_ij.
      
      int k; // for iterating over 6 constraints for a given link
      for (k=0; k<6; k++) {
	// part k of piecewise linear constraint: phi_ij >= f(y_ij)
	aptrb[row] = subindex;
	for (c=0; c<C; c++) {
	  asub[subindex] = i*J*C + j*C + c; // which variable we're indexing. Here, it's beta_ijc.
	  aval[subindex] = 0.0;
	  for (n=0; n<N; n++)
	    aval[subindex] += (y_coeff[k] * p.vol[c*N + n] * alpha[i*C*N + c*N + n]);
	  subindex++;
	} // for each client utilizing the link
	
	asub[subindex] = I*J*C + i*J + j;
	aval[subindex] = -1.0; // for -phi_ij in the constraint
	subindex ++;

	// finished constraint for k'th part
	aptre[row] = subindex;
	row ++; 
      } // end of 6 constraints per link
    } 
  } // end of for each link

  // traffic inclusion constraints
  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) { // for each constraint
      aptrb[row] = subindex;
      for (j=0; j<J; j++) {
	asub[subindex] = i*J*C + j*C + c;
	aval[subindex] = 1;
	subindex ++;
      }
      aptre[row] = subindex;
      row ++;
    }
  }

  //printf("Routing: Number of constraints %d %d number of nonzero coeffs %d %d\n", row, 6*I*J+I*C, subindex, 7*I*J*C + 6*I*J);

  return 0;

}

int prepare_constbounds_dc(struct params_t p, double* alpha, MSKboundkeye* bkc, double* blc, double* buc) {
  
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  int row = 0;

  int i, j, c, n;

  // load balancing constraint bounds - 6 bounds per link
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) { // for every link      
      int k;
      for (k=0; k<6; k++) {
	// part k of piecewise linear function
	bkc[row] = MSK_BK_UP;
	blc[row] = -MSK_INFINITY;
	buc[row] = cap_coeff[k] * p.cap[i*J + j];
	row ++;
      } // for each link constraint
    }
  } // for each link

  // traffic inclusion constraints bounds
  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) {
      bkc[row] = MSK_BK_FX;
      blc[row] = 1.0;
      buc[row] = 1.0;
      row ++;
    }
  }

  return 0;

}

int prepare_varbounds_dc(struct params_t p, MSKboundkeye* bkx, double* blx, double* bux) {
  
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  int* lnum = p.lnum;

  int i, j, c, n;

  int var = 0;

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) { // for every beta_ijc variable
	if (j < lnum[i]) {	  
	  if (p.perf[i*J*C + j*C + c] >= 0.000001) { // valid performance information on this link
 	    bkx[var] = MSK_BK_LO;
	    blx[var] = 0.0;
	    bux[var] = MSK_INFINITY;
	  }
	  else { // assume not reachable, can't use this link
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
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) { // for every phi_ij variable
      bkx[var] = MSK_BK_LO;
      blx[var] = 0.0;
      bux[var] = MSK_INFINITY;
      var ++;
    }
  }

  return 0;

}

int setup_dc(struct params_t p, struct msk_problem_t* dc, double* alpha) {
  
  // prepare problem related arrays
  // 1. objective coefficients
  // 2. constraint matrix sparse description arrays
  // 3. constraint bound arrays
  // 4. update provided structure with problem information

  prepare_objective_dc(p, alpha, dc->obj_coeff);
  prepare_constraints_dc(p, alpha, dc->aptrb, dc->aptre, dc->asub, dc->aval);
  prepare_constbounds_dc(p, alpha, dc->bkc, dc->blc, dc->buc);
  prepare_varbounds_dc(p, dc->bkx, dc->blx, dc->bux);

  return 0;

}

// allocate msk_problem_t structure for DC "master" problem
struct msk_problem_t* alloc_dc(struct params_t p) {
  
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  // coefficients of the objective function
  // -- corresponds to number of variables involved
  double* obj_coeff = (double*)malloc((I*J*C + I*J)*sizeof(double));
  
  // constraints
  int* aptrb = (int*)malloc((6*I*J+I*C)*sizeof(int)); // number of rows
  int* aptre = (int*)malloc((6*I*J+I*C)*sizeof(int)); // number of rows
  int* asub  = (int*)malloc((7*I*J*C + 6*I*J) * sizeof(int)); // number of nonzero terms
  double* aval = (double*)malloc((7*I*J*C + 6*I*J) * sizeof(double)); // number of nonzero terms

  // bounds on constraints - depends on number of constraints
  MSKboundkeye* bkc = (MSKboundkeye*)malloc((6*I*J+I*C)*sizeof(MSKboundkeye));
  double* blc = (double*)malloc((6*I*J+I*C)*sizeof(double));
  double* buc = (double*)malloc((6*I*J+I*C)*sizeof(double));

  // bounds on variables
  MSKboundkeye* bkx = (MSKboundkeye*)malloc((I*J*C + I*J)*sizeof(MSKboundkeye));
  double* blx = (double*)malloc((I*J*C + I*J)*sizeof(double));
  double* bux = (double*)malloc((I*J*C + I*J)*sizeof(double));

  // Memory alloc checking
  if(obj_coeff == 0 || aptrb == 0 || aptre == 0 || asub == 0 || aval == 0 || bkc == 0 || blc == 0 || buc == 0 || 
       bkx == 0 || blx == 0 || bux == 0) {
    printf("Memory allocation error in data center subproblem!\n");
    exit(0);
  }

  struct msk_problem_t *mp = (struct msk_problem_t*)malloc(sizeof(struct msk_problem_t));
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
int dealloc_dc(struct msk_problem_t* dc) {

  free(dc->obj_coeff);
  free(dc->aptrb);
  free(dc->aptre);
  free(dc->asub);
  free(dc->aval);
  free(dc->bkc);
  free(dc->blc);
  free(dc->buc);
  free(dc->bkx);
  free(dc->blx);
  free(dc->bux);
  free(dc);
  
}

// prints out A vectors for inspection for huge a_ij problem
int check_avec(int* asub, double* aval, int size) {
  int i;
  for (i=0; i<size; i++) {
    printf("%d %lf\n", asub[i], aval[i]);
  }
  return 0;
}

int optimize_dc(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* beta, MSKrealt* objvalue, int print_sol) {

  MSKtask_t    task = NULL;
  MSKidxt i,j;
  MSKrescodee  r;
  
  int numvar = p.I * p.J * p.C + p.I * p.J;
  int numcon = 6 * p.I * p.J + p.I * p.C;
  int numanz = 7 * p.I * p.J * p.C + 6 * p.I * p.J;

  //check_avec(mp->asub, mp->aval, numanz);

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
			     MSK_SOL_BAS,
			     NULL,
			     &solsta);
      switch(solsta) {
      case MSK_SOL_STA_OPTIMAL:   
      case MSK_SOL_STA_NEAR_OPTIMAL:

	// get complete primal solution
	MSK_getsolutionslice(task,
			     MSK_SOL_BAS,    /* Request the basic solution. */
			     MSK_SOL_ITEM_XX,/* Which part of solution.     */
			     0,              /* Index of first variable.    */
			     numvar,         /* Index of last variable+1.   */
			     beta);
	
	if (print_sol > 0) {
	  printf("Optimal primal solution\n");
	  /* for(j=0; j<numvar; ++j) */
	  /*   printf("beta[%d]: %e\n", j, beta[j]); */
	}

	MSK_getprimalobj(task, MSK_SOL_BAS, objvalue);
	
	break;
      case MSK_SOL_STA_DUAL_INFEAS_CER:
      case MSK_SOL_STA_PRIM_INFEAS_CER:
      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
	printf("Routing problem: Primal or dual infeasibility certificate found.\n");
	// write problem data to a file for future inspection
	MSK_writedata(task, "problem_beta.opf");
	MSK_deletetask(&task);
	rvalue = 1;
	break;
        
      case MSK_SOL_STA_UNKNOWN:
	printf("Routing problem: The status of the solution could not be determined.\n");
	// write problem data to a file for future inspection
	MSK_writedata(task, "problem_beta.opf");
	rvalue = 1;
	break;
      default:
	printf("Routing problem: Other solution status.");
	// write problem data to a file for future inspection
	MSK_writedata(task, "problem_beta.opf");
	rvalue = 1;
	break;
      }
    }
    else {
      printf("Routing problem: Error while optimizing.\n");
      // write problem data to a file for future inspection
      MSK_writedata(task, "problem_beta.opf");
      rvalue = 1;
    }
  }
  
  if (r != MSK_RES_OK) {
    /* In case of an error print error code and description. */      
    char symname[MSK_MAX_STR_LEN];
    char desc[MSK_MAX_STR_LEN];
    
    printf("Routing problem: An error occurred while optimizing.\n");
    MSK_getcodedesc (r,
		     symname,
		     desc);
    printf("Error %s - '%s'\n",symname,desc);
    
    rvalue = 1;

  }
  
  MSK_deletetask(&task);

  return rvalue;

}

// initialize beta such that one random link per data center is picked
int init_beta(struct params_t p, double* beta) {

  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  
  int i, j, c, n;

  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) {

      // strategy 1. pick a random outgoing link to start with

      // for each DC and each client
      // initially set all betas to be zero
      for (j=0; j<J; j++) {
      	beta[i*J*C + j*C + c] = 0.0;
      }
      // randomly pick a link to send traffic out of
      // for this client
      j = (int)((double)rand() * J/RAND_MAX);
      beta[i*J*C + j*C + c] = 1.0;

      // strategy 2. equitable distribution

      /* for (j=0; j<J; j++) */
      /* 	beta[i*J*C + j*C + c] = 1.0/J; */
    }
  }
  
  return 0;

}
