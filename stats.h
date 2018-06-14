#pragma once

struct msk_stats_t {
  int *alpha_final_cdf, *beta_final_cdf, *alpha_interim_cdf, *beta_interim_cdf, *numiter_cdf;
  double *objdiff_values;
  int alpha_size, beta_size, objdiff_runs;
};

// function declarations
int write_cdfs(struct msk_stats_t *s);
int alloc_cdf(struct msk_stats_t* s, int I, int J, int C, int N, int runs);
int dealloc_cdf(struct msk_stats_t* s);
int updatecdf_finals(struct msk_stats_t *s, double* alpha, double* beta, int numiter);
int updatecdf_interim(struct msk_stats_t *s, double* alpha, double* beta);
