#pragma once

#include "mosek.h"

// multiplicative factor on link penalty to get it to 250ms
#define SC_PHI 23.4375

// fraction of capacity used to bound link congestion cost
#define CAPFRAC 0.95
#define CAPFRAC2 0.945

#define GLOBAL_PROBLEM_FILENAME "foo.opf"

// maximum number of primal decomposition iterations
#define MAXITER 200

double total_traffic(double* vol, int size);
int init_mosek_env(MSKenv_t* env);

int y_coeff[6];
double cap_coeff[6];
