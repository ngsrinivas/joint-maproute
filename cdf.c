#include "cdf.h"
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

// Decided to implement a new CDF module, because a lot of things
// are much more general here, as opposed to stats.c, which only
// has a reasonably general int cdf implementation.

struct cdf_t* cdf_init_doublecdf(int totalitems, int numbins) {

  struct cdf_t* new_cdf = (struct cdf_t*)malloc(sizeof(struct cdf_t));
  if (new_cdf == NULL) { fprintf(stderr, "CDF: Could not initialize new_cdf!\n"); exit(0); }

  new_cdf->totalitems = totalitems;
  new_cdf->items = (double*)malloc(totalitems * sizeof(double));
  if (new_cdf->items == NULL) { fprintf(stderr, "CDF: Could not initialize item list! %d\n", totalitems); exit(0); }
  
  new_cdf->bins = (double*)calloc(numbins + 1, sizeof(double));
  if (new_cdf->bins == NULL) { fprintf(stderr, "CDF: Could not initialize bins %d\n", numbins); exit(0); }
  
  new_cdf->numitems = 0;
  new_cdf->numbins = numbins;
  new_cdf->minvalue = new_cdf->maxvalue = 0.0;
  return new_cdf;
}

// Tear down a CDF structure
int cdf_dealloc_doublecdf(struct cdf_t* cdf) {
  free(cdf->items);
  free(cdf->bins);
  free(cdf);
  return 0;
}

// Add an item to be included in the CDF list
int cdf_additem(struct cdf_t* cdf, double item) {
  if (cdf->numitems < cdf->totalitems) { // valid operation
    cdf->items[cdf->numitems] = item;
    cdf->numitems ++;
    // update min and max values
    if (cdf->numitems == 1) {
      cdf->minvalue = item;
      cdf->maxvalue = item;
    } else {
      if (cdf->minvalue > item) cdf->minvalue = item;
      if (cdf->maxvalue < item) cdf->maxvalue = item;
    }
  }
  else {
    fprintf(stderr, "CDF: Item list not big enough! %d %lf\n", cdf->totalitems, item);
    exit(0);
  }
  return 0;
}

// compute CDF and store values in bins.
int cdf_compute_doublecdf(struct cdf_t* cdf) {
  if (cdf->numitems > 1) {
    double range = cdf->maxvalue - cdf->minvalue;
    if (range >= 0.000001) {
      // bin range into numbins
      int i, binindex;
      double multfactor = cdf->numbins * 1.0/range;
      for (i=0; i<cdf->numitems; i++) {
	binindex = (int)(ceil((cdf->items[i] - cdf->minvalue) * multfactor));
	// fix for one-off floating point problems
	if (binindex == cdf->numbins + 1) binindex = cdf->numbins;
	assert(binindex >= 0 && binindex <= cdf->numbins);
	cdf->bins[binindex] += 1;
      }

      // at this stage, this is still a frequency count of the various ranges
      double cdf_value = 0.0;
      for (i=0; i<=cdf->numbins; i++) {
	cdf_value += cdf->bins[i]; // cumulative frequency
	cdf->bins[i] = cdf_value/cdf->numitems;
      }
      // At this point, cdf->bins[i] contains the CDF, ie, P(X <= x) for every x value, where the xvalue for bins[i] is minvalue + i/multfactor.
      return 0;
    }
    else {
      fprintf(stderr, "CDF: Range of item list not big enough to bin for CDF %lf\n", range);
      return 1;
    }
  }
  else {
    fprintf(stderr, "CDF: Item list not big enough to compute CDF\n");
    return 1;
  }
}

// To be called with an open writeabe file descriptor AFTER calling
// compute_doublecdf on the cdf_t object.
int cdf_writecdf(FILE* fp, struct cdf_t* cdf) {

  double multfactor = (cdf->maxvalue - cdf->minvalue)/cdf->numbins;
  
  int i;
  for (i=0; i<=cdf->numbins; i++) {
    // x, P(X <= x)
    fprintf(fp, "%lf %lf\n", cdf->minvalue + (i*multfactor), cdf->bins[i]);
  }
  
  return 0;
}

