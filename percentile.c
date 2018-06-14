#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define TOL 0.00001

int main(int argc, char** argv) {
  
  if (argc < 3) {
    fprintf(stderr, "Usage: cat cdf-datafile | %s <variance-value> percentile1 [percentile2 ...]\n", argv[0]);
    fprintf(stderr, "Percentile input values must be sorted in increasing order for correctness. Percentiles must be integers.\n");
    exit(0);
  }
  
  double reqvar = atof(argv[1]);
  double npercs = argc - 2;
  int* percvalues = (int*)malloc(npercs * sizeof(int));

  int i;
  // Read in percentile numbers, and print a header for the 
  // output along the way
  for (i=0; i < npercs; i++) {
    percvalues[i] = atoi(argv[i+2]);
  }

  double x, p, prev_x = 0.0, prev_p = 0.0;
  int lastperc_index = 0;
  double required_percentile = percvalues[0] / 100.0;
  
  printf("%lf ", reqvar);

  while (lastperc_index < npercs) {
    if (scanf("%lf %lf\n", &x, &p) > 0) { // just read a line of data
      while (lastperc_index < npercs && p >= required_percentile) {
	printf("%lf ", x);	
	lastperc_index ++;
	if (lastperc_index < npercs) required_percentile = percvalues[lastperc_index] / 100.0;
      }
    }
    else break;
  }
  
  printf("\n");
  
  /* // old version of the loop */
  /* while (scanf("%lf %lf\n", &x, &p) > 0) { */
  /*   if (p >= required_percentile) { */
  /*     // print out data point */
  /*     printf("%lf %lf\n", reqvar, x); */
      
  /*     lastperc_index ++; */
  /*     if (lastperc_index < npercs) { */
  /* 	required_percentile = percvalues[lastperc_index]; */
  /*     } */
  /*     else { */
  /* 	break; */
  /*     } */
  /*   } */
  /* } */

  return 0;

}

