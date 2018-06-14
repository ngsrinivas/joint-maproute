#pragma once

#include <stdio.h>
#include <stdlib.h>

struct cdf_t {
  int numitems, totalitems, numbins;
  double* items;
  double* bins;
  double minvalue, maxvalue;
};

struct cdf_t* cdf_init_doublecdf(int totalitems, int numbins);
int cdf_dealloc_doublecdf(struct cdf_t* cdf);
int cdf_additem(struct cdf_t* cdf, double item);
int cdf_writecdf(FILE *fp, struct cdf_t* cdf);
int cdf_compute_doublecdf(struct cdf_t* cdf);
