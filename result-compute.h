#pragma once

#include "params.h"

// X based functions to print various numbers
// of the optimization problem.
double compute_perf_X(struct params_t* p, double* X);
double compute_total_price_X(struct params_t* p, double* X);
double compute_price_per_byte_X(struct params_t* p, double* X);
double compute_just_wrtt_X(struct params_t* p, double* X);

int print_link_utilizations_X(struct params_t* p, double* X);

// similar alpha beta based functions.
double compute_just_wrtt_alphabeta(struct params_t* p, double* alpha, double* beta);
double compute_perf_alphabeta(struct params_t* p, double* alpha, double* beta);
double compute_perf2_alphabeta(struct params_t* p, double* alpha, double* beta);
double compute_total_price_alphabeta(struct params_t* p, double* alpha, double* beta);
double compute_price_per_byte_alphabeta(struct params_t* p, double* alpha, double* beta);
double compute_objperf_X(struct params_t* p, double* X);
double compute_objperf_alphabeta(struct params_t* p, double* alpha, double* beta);

double compute_wrtt_with_backbone(struct params_t* p, double* alpha, double* X, double* blats);

// per client computation functions
double compute_avgcost_c(struct params_t p, double* alpha, double* beta, int c);
double compute_avglats_c(struct params_t p, double* alpha, double* beta, int c);

