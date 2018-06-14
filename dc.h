#pragma once
#include "params.h"

int setup_dc(struct params_t p, struct msk_problem_t* dc, double* alpha);
struct msk_problem_t* alloc_dc(struct params_t p);
int dealloc_dc(struct msk_problem_t* dc);
int optimize_dc(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* beta, MSKrealt* objvalue, int print_sol);
int init_beta(struct params_t p, double* beta);

int setup_dc_lat(struct params_t p, struct msk_problem_t* dc, double* alpha);
struct msk_problem_t* alloc_dc_lat(struct params_t p);
int dealloc_dc_lat(struct msk_problem_t* dc);
int optimize_dc_lat(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* beta, MSKrealt* objvalue, int print_sol);

int setup_dc_95_lat(struct params_t p, struct msk_problem_t* dc, double* alpha);
struct msk_problem_t* alloc_dc_95_lat(struct params_t p);
int dealloc_dc_95_lat(struct msk_problem_t* dc);
int optimize_dc_95_lat(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* beta, MSKrealt* objvalue, int print_sol);

int setup_dc_cost(struct params_t p, struct msk_problem_t* dc, double* alpha);
struct msk_problem_t* alloc_dc_cost(struct params_t p);
int dealloc_dc_cost(struct msk_problem_t* dc);
int optimize_dc_cost(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* beta, MSKrealt* objvalue, int print_sol);

