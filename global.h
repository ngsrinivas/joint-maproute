#pragma once

struct msk_problem_t* alloc_global(struct params_t p);
int setup_global(struct params_t p, struct msk_problem_t *mp);
int optimize_global(MSKenv_t* env, struct params_t p, struct msk_problem_t* mp, double* X, MSKrealt* objvalue, int print_sol);


