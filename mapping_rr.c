#include "helper.h"
#include<stdlib.h>
#include "params.h"
#include "mosek.h"

/*
  
  Round robin mapping. Doesn't care about capacities or 
  the current routing.
  
 */


int optimize_mapping_rr(struct params_t p, double* alpha) {
  
  int i, j, c, n, I, J, C, N;
  I = p.I;
  J = p.J;
  C = p.C;
  N = p.N;

  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) {
      for (n=0; n<N; n++) {
	alpha[i*C*N + c*N + n] = 1.0/I;
      }
    }
  }
  
  return 0;

}
