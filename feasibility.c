#include<stdio.h>
#include "params.h"

int check_feasibility(struct params_t p, double* alpha, double* beta) {
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;
  double K = p.K;
  double mu = p.mu;
  double TT = p.TT;

  int i, j, c, n;

  int success = 1;
  double tolerance = 0.00001;
  
  // set 1. load balancing constraints
  for (i=0; i<I; i++) { // for each constraint
    double load_i = 0.0;
    for (c=0; c<C; c++) {
      for (n=0; n<N; n++) {
	load_i += (alpha[i*C*N+c*N+n] * (p.vol[c*N+n]/TT));
      }
    }

    // verify the current inequalities directly..
    double load_diff = load_i - p.w[i];

    if (load_i > p.w[i] + p.eps[i] + tolerance) {
      fprintf(stderr, "***Set 1 constraint %d upper bound violated %lf\n", i, load_i - p.w[i] - p.eps[i]);
      success = 0;
    }
    if (load_i < p.w[i] - p.eps[i] - tolerance) {
      fprintf(stderr, "***Set 1 constraint %d lower bound violated %lf\n", i, p.w[i] - p.eps[i] - load_i);
      success = 0;
    }
  }

  // set 2. link capacity constraints
  // Commented out, because they are no more constraints.
  /* for (i=0; i<I; i++) { */
  /*   for (j=0; j<J; j++) { // for each constraint */
  /*     double load_ij = 0.0; */
  /*     for (c=0; c<C; c++) { */
  /* 	for (n=0; n<N; n++) { */
  /* 	  load_ij += (alpha[i*C*N+c*N+n] * beta[i*J*C+j*C+c] * p.vol[c*N+n]); */
  /* 	} */
  /*     } */
  /*     if (load_ij > mu * p.cap[i*J+j] + tolerance) { */
  /* 	fprintf(stderr, "***Set 2 constraint %d violated %lf\n", i*j, load_ij - (mu * p.cap[i*J+j])); */
  /* 	success = 0; */
  /*     } */
  /*   } */
  /* } */

  // set 3. traffic inclusion constraints
  for (c=0; c<C; c++) {
    for (n=0; n<N; n++) { // for each constraint
      double alpha_total = 0.0;
      for (i=0; i<I; i++)
	alpha_total += alpha[i*C*N+c*N+n];
      if((alpha_total - 1.00 >= tolerance) || (1.00 - alpha_total >= tolerance)) {
	fprintf(stderr, "***Set 3a constraint %d violated %lf\n", c*N+n, alpha_total);
	success = 0;
      }
    }
  }
  
  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) { // for each constraint
      double beta_total = 0.0;
      for (j=0; j<J; j++) 
	beta_total += beta[i*J*C+j*C+c];
      if((beta_total - 1.00 >= tolerance) || (1.00 - beta_total >= tolerance)) {
	fprintf(stderr, "***Set 3b constraint %d violated %lf\n", i*C+c, beta_total);
	success = 0;
      }
    }
  }

  if (success) {
    //fprintf(stderr, "####################################\n");
    //fprintf(stderr, "### All feasibility checks succeeded\n");
    //fprintf(stderr, "####################################\n");
  }

  return success; // 1 if feasible, 0 if not
  
}
