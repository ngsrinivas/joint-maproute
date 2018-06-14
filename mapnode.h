#pragma once

#include "params.h"
#include "mosek.h"

int optimize_mapnode(MSKenv_t env, struct params_t p, struct msk_problem_t* mp, double* alpha, MSKrealt* objvalue, int print_sol);
int dealloc_mapnode(struct msk_problem_t *mp);
struct msk_problem_t* alloc_mapnode(struct params_t p);
int setup_mapnode(struct params_t p, struct msk_problem_t *mp, double* beta);

int optimize_mapnode_rtaware(MSKenv_t env, struct params_t p, struct msk_problem_t* mp, double* alpha, MSKrealt* objvalue, int print_sol);
int dealloc_mapnode_rtaware(struct msk_problem_t *mp);
struct msk_problem_t* alloc_mapnode_rtaware(struct params_t p);
int setup_mapnode_rtaware(struct params_t p, struct msk_problem_t *mp, double* beta);

int optimize_mapnode_95_rtaware(MSKenv_t env, struct params_t p, struct msk_problem_t* mp, double* alpha, MSKrealt* objvalue, int print_sol);
int dealloc_mapnode_95_rtaware(struct msk_problem_t *mp);
struct msk_problem_t* alloc_mapnode_95_rtaware(struct params_t p);
int setup_mapnode_95_rtaware(struct params_t p, struct msk_problem_t *mp, double* beta);

int optimize_mapnode_capped_closest(MSKenv_t env, struct params_t p, struct msk_problem_t* mp, double* alpha, MSKrealt* objvalue, int print_sol);
int dealloc_mapnode_capped_closest(struct msk_problem_t *mp);
struct msk_problem_t* alloc_mapnode_capped_closest(struct params_t p);
int setup_mapnode_capped_closest(struct params_t p, struct msk_problem_t *mp);

int optimize_mapnode_95_capped_closest(MSKenv_t env, struct params_t p, struct msk_problem_t* mp, double* alpha, MSKrealt* objvalue, int print_sol);
int dealloc_mapnode_95_capped_closest(struct msk_problem_t *mp);
struct msk_problem_t* alloc_mapnode_95_capped_closest(struct params_t p);
int setup_mapnode_95_capped_closest(struct params_t p, struct msk_problem_t *mp);

int optimize_mapping_rr(struct params_t p, double* alpha);

int optimize_mapping_closest(struct params_t p, double* alpha);

int init_alpha(struct params_t p, double* alpha);



