#include "params.h"
#include "helper.h"
#include <assert.h>

double compute_just_wrtt_X(struct params_t* p, double* X) {
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	double vol_c = 0.0;
	for (n=0; n<N; n++)
	  vol_c += p->vol[c*N + n];
	value += (vol_c/TT * p->perf[i*J*C + j*C + c] * X[i*J*C + j*C + c]);	
      }
    }
  }
  
  return value;
}

double compute_objperf_X(struct params_t* p, double* X) {
  // return performance per byte, including phi_ij variables.
  
  int I = p->I;
  int J = p->J;
  int C = p->C;

  int i, j;
  double value = 0.0;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p->cap[i*J + j] >= 0.000001) {
	value += (X[I*J*C + i*J + j]/(p->cap[i*J + j] * p->total_links));
      }
    }
  }
  
  value *= SC_PHI;
  value += compute_just_wrtt_X(p, X);
  
  return value;
  
}

double compute_perf_X(struct params_t* p, double* X) {
  // return performance per byte, including phi_ij variables.
  
  fprintf(stderr, "compute_perf_X assumes existence of phi_ij variables in the provided solution X\n");

  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p->cap[i*J+j] >= 0.000001) {
	for (c=0; c<C; c++) {
	  double vol_c = 0.0;
	  for (n=0; n<N; n++)
	    vol_c += p->vol[c*N + n];
	  value += ((vol_c/p->TT) * X[i*J*C + j*C + c] * (SC_PHI * X[I*J*C + i*J + j]/(p->cap[i*J + j])));
	}
      }
    }
  }
  
  value += compute_just_wrtt_X(p, X);
  
  return value;
}

double compute_total_price_X(struct params_t* p, double* X) {
  // return total cost, including phi_ij variables.

  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	double vol_c = 0.0;
	for (n=0; n<N; n++)
	  vol_c += p->vol[c*N + n];
	value += (vol_c * p->price[i*J + j] * X[i*J*C + j*C + c]);
      }
    }
  }
  
  return value;
}

double compute_price_per_byte_X(struct params_t* p, double* X) {
  return compute_total_price_X(p, X)/p->TT;
}

int print_link_utilizations_X(struct params_t* p, double* X) {
  
  int i, j, c, n, I, J, C, N;
  I = p->I;
  J = p->J;
  C = p->C;
  N = p->N;
  
  for (i=0; i<I; i++) {
    for (j=0; j<p->lnum[i]; j++) {
      double load = 0.0;
      for (c=0; c<C; c++) {
	double vol_c = 0.0;
	for (n=0; n<N; n++)
	  vol_c += p->vol[c*N + n];
	load += (vol_c * X[i*J*C + j*C + c]);
      }
      printf("%d %d %lf\n", i, j, load/p->cap[i*J+j]);
    } // for every link
  }

  return 0;

}

double compute_just_wrtt_alphabeta(struct params_t* p, double* alpha, double* beta) {
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	for (n=0; n<N; n++) {
	  value += (p->vol[c*N + n]/TT * p->perf[i*J*C + j*C + c] * beta[i*J*C + j*C + c] * alpha[i*C*N + c*N + n]);
	}
      }
    }
  }
  
  return value;
}

double compute_phi(double y, double cap, double scaling) {
  
  assert(cap >= 0.000001);
  double frac = y/cap;
  double out;
  
  if (frac < (1.0/3.0)) 
    out = y;
  else if (frac < (2.0/3.0)) 
    out = 3*y - (2.0/3.0)*cap;
  else if (frac < 0.9)
    out = 10*y - (16.0/3.0)*cap;
  else if (frac < 1)
    out = 70*y - (178.0/3.0)*cap;
  else out = 5000*y - (14968.0/3.0)*cap;
  
  return out * scaling / cap;

}

double compute_perf2_alphabeta(struct params_t* p, double* alpha, double* beta) {
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;

  // Implementation for schemes where routing clearly has the phi_ij in it
  // otherwise, beta[IJC + iJ + j] would not exist!
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p->cap[i*J+j] >= 0.000001) {
  	for (c=0; c<C; c++) {
  	  for (n=0; n<N; n++) {
  	    value += ((p->vol[c*N + n]/p->TT) * beta[i*J*C + j*C + c] * alpha[i*C*N + c*N + n] * (SC_PHI * beta[I*J*C + i*J + j]/(p->cap[i*J + j])));
  	  }
  	}
      }
    }
  }
  
  value += compute_just_wrtt_alphabeta(p, alpha, beta);
  
  return value;
}

double compute_perf_alphabeta(struct params_t* p, double* alpha, double* beta) {
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;

  // recomputing phi_ij just to be safe.. and have to do it here
  // anyway because some routing schemes don't have these variables now.
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      // determine if legitimate link
      if (p->cap[i*J + j] >= 0.000001) {
	double vol_cn = 0.0;
	for (c=0; c<C; c++) {
	  for (n=0; n<N; n++) {
	    vol_cn += (p->vol[c*N + n] * alpha[i*C*N + c*N + n] * beta[i*J*C + j*C + c]);
	  }
	}	
	value += ((vol_cn/p->TT) * compute_phi(vol_cn, p->cap[i*J + j], SC_PHI));
      } // valid link
    }
  } // every link
	
  // Implementation for schemes where routing clearly has the phi_ij in it
  // otherwise, beta[IJC + iJ + j] would not exist!
  /* for (i=0; i<I; i++) { */
  /*   for (j=0; j<J; j++) { */
  /*     if (p->cap[i*J+j] >= 0.000001) { */
  /* 	for (c=0; c<C; c++) { */
  /* 	  for (n=0; n<N; n++) { */
  /* 	    value += ((p->vol[c*N + n]/p->TT) * beta[i*J*C + j*C + c] * alpha[i*C*N + c*N + n] * (SC_PHI * beta[I*J*C + i*J + j]/(p->cap[i*J + j]))); */
  /* 	  } */
  /* 	} */
  /*     } */
  /*   } */
  /* } */
  
  value += compute_just_wrtt_alphabeta(p, alpha, beta);
  
  return value;
}

double compute_objperf_alphabeta(struct params_t* p, double* alpha, double* beta) {
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;

  // recomputing phi_ij just to be safe.. and have to do it here
  // anyway because some routing schemes don't have these variables now.
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      // determine if legitimate link
      if (p->cap[i*J + j] >= 0.000001) {
	double vol_cn = 0.0;
	for (c=0; c<C; c++) {
	  for (n=0; n<N; n++) {
	    vol_cn += (p->vol[c*N + n] * alpha[i*C*N + c*N + n] * beta[i*J*C + j*C + c]);
	  }
	}	
	value += (compute_phi(vol_cn, p->cap[i*J + j], SC_PHI) / p->total_links);
      } // valid link
    }
  } // every link
	
  
  value += compute_just_wrtt_alphabeta(p, alpha, beta);
  
  return value;
}


double compute_total_price_alphabeta(struct params_t* p, double* alpha, double* beta) {
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  double L = p->L;
  double K = p->K;
  double TT = p->TT;
  
  int i, j, c, n;
  
  double value = 0.0;
  
  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      for (c=0; c<C; c++) {
	double vol_c = 0.0;
	for (n=0; n<N; n++) {
	  value += (p->vol[c*N+n] * p->price[i*J + j] * beta[i*J*C + j*C + c] * alpha[i*C*N + c*N + n]);
	}
      }
    }
  }
  
  return value;

}

double compute_price_per_byte_alphabeta(struct params_t* p, double* alpha, double* beta) {
  return compute_total_price_alphabeta(p, alpha, beta)/p->TT;
}

// Compute "actual" wRTT for entact, given alpha, X, backbone latencies
// and some basic problem information (params_t).
double compute_wrtt_with_backbone(struct params_t* p, double* alpha, double* X, double* blats) {
  
  int i, j, c, n, e;
  int I = p->I;
  int J = p->J;
  int C = p->C;
  int N = p->N;
  
  // to store the best performant RTT from i to c.
  double* perf = (double*)malloc(I*C*sizeof(double));
  
  // compute perf_ic, the best round trip latencies for 
  // paths from i to c.
  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) {
      double minlat = 100000;
      for (j=0; j<J; j++) {
	double curr_perf = p->perf[i*J*C + j*C + c];
	if (curr_perf >= 0.000001 && curr_perf < minlat) minlat = curr_perf;
      }
      perf[i*C + c] = minlat;
    }
  }

  // compute wRTT for this set of paths, using backbone, perf_ic
  // and perf_ijc values.
  double value = 0.0;
  for (i=0; i<I; i++) {
    for (e=0; e<I; e++) {
      for (j=0; j<J; j++) {
	for (c=0; c<C; c++) {
	  for (n=0; n<N; n++) {
	    value += (alpha[i*C*N + c*N + n] * (p->vol[c*N + n]/p->TT) * X[e*J*C + j*C + c] * (blats[i*I + e] + (perf[i*C + c]/2.0) + (p->perf[i*J*C + j*C + c]/2.0)));
	  }
	}
      }
    }
  }
  
  return value;

}

double compute_avglats_c(struct params_t p, double* alpha, double* beta, int c) {
  
  int i, j, n;
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  double avglats = 0.0;
  double vol_c = 0.0;
  
  for (n=0; n<N; n++) vol_c += p.vol[c*N + n];
  
  //printf("vol_c for client %d is %lf\n", c, vol_c);

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p.cap[i*J + j] >= 0.000001 && p.perf[i*J*C + j*C + c] >= 0.000001) {
	for (n=0; n<N; n++) {
	  avglats += (alpha[i*C*N + c*N + n] * (p.vol[c*N + n]/vol_c) * beta[i*J*C + j*C + c] * (p.perf[i*J*C + j*C + c] + (SC_PHI * beta[I*J*C + i*J + j]/p.cap[i*J + j])));
	  //printf("alpha_icn %lf beta_ijc %lf phi_ij %lf perf_ijc %lf vol_cn %lf cap_ij %lf\n", alpha[i*C*N + c*N + n], beta[i*J*C + j*C + c], beta[I*J*C + i*J + j], p.perf[i*J*C + j*C + c], p.vol[c*N + n], p.cap[i*J + j]);
	}
      } // valid link
    }
  }
  
  if (avglats != avglats) {// NaN test
    printf("NaN trouble! client number %d\n", c);
    exit(0);
  }

  return avglats;

}

double compute_avgcost_c(struct params_t p, double* alpha, double* beta, int c) {
  
  int i, j, n;
  int I = p.I;
  int J = p.J;
  int C = p.C;
  int N = p.N;

  double avgcost = 0.0;
  double vol_c = 0.0;
  
  for (n=0; n<N; n++) vol_c += p.vol[c*N + n];

  for (i=0; i<I; i++) {
    for (j=0; j<J; j++) {
      if (p.cap[i*J + j] >= 0.000001 && p.perf[i*J*C + j*C + c] >= 0.000001) {
	for (n=0; n<N; n++) {
	  avgcost += (alpha[i*C*N + c*N + n] * (p.vol[c*N + n]/vol_c) * beta[i*J*C + j*C + c] * p.price[i*J + j]);
	}
      } // valid link
    }
  }

  double gbs_per_req = 65.8 * 0.000001;
  return avgcost / gbs_per_req;

}


