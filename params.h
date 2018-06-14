#pragma once
#include "mosek.h"

struct params_t {
  double *price, *perf, *vol, *w, *eps, *cap;
  double mu, K, TT, L;
  int I, J, C, N, total_links;
  int* lnum; // link numbers per DC
};

struct msk_problem_t {
  double *obj_coeff, *aval, *blc, *buc, *blx, *bux;
  int *aptrb, *aptre, *asub;
  MSKboundkeye *bkx, *bkc;
};

struct msk_qconst_t {
  int *qsubk, *qsubi, *qsubj;
  double *qval;
};


