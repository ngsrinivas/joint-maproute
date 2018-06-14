#include "helper.h"
#include<stdlib.h>
#include "params.h"
#include "mosek.h"

/*
  
  Mapping performed by routing to the DC which gives the 
  path with the smallest latency. But don't care about 
  DC capacity or link capacities there, or the current routing.
  
 */

int optimize_mapping_closest(struct params_t p, double* alpha) {
  
  int i, j, c, n, I, J, C, N;
  I = p.I;
  J = p.J;
  C = p.C;
  N = p.N;
  
  // set all alpha to zero to start off with
  for (i=0; i<I; i++) {
    for (c=0; c<C; c++) {
      for (n=0; n<N; n++) {
	alpha[i*C*N + c*N + n] = 0.0;
      }
    }
  }

  for (c=0; c<C; c++) {
    
    double minlat = MSK_INFINITY;
    int min_i = -1;

    // check best performing link
    for (i=0; i<I; i++) {
      for (j=0; j<J; j++) {
	double curr_perf = p.perf[i*J*C + j*C + c];
	if (curr_perf >= 0.000001 && curr_perf < minlat) {
	  minlat = p.perf[i*J*C + j*C + c];
	  min_i = i;
	}
      }
    }
    
    for (n=0; n<N; n++)
      alpha[min_i*C*N + c*N + n] = 1.0;
  } // for every client
  
  return 0; // no problem!

}

