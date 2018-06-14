#include<stdlib.h>
#include<stdio.h>

#define MAXINDEX 10

FILE* fps[MAXINDEX];
int last_index = 0;

int diag_init(char* filename) {
  if (last_index >= MAXINDEX) {
    printf("Diagnostics needs more file pointer space.\n");
    exit(0);
  }
  else {
    fps[last_index] = fopen(filename, "w");
    if (fps != 0) {
      last_index ++;
      return last_index-1;
    }
    else {
      printf("Diagnostics could not open new file for writing!\n");
      exit(0);
    }
  }
}

int diag_write(int fpindex, double x, double y) {
  fprintf(fps[fpindex], "%lf %lf\n", x, y);
  return 0;
}

int diag_write_list(int fpindex, double x, double* y, int size) {
  fprintf(fps[fpindex], "%lf", x);
  int i;
  for (i=0; i<size; i++)
    fprintf(fps[fpindex], " %lf", y[i]);
  fprintf(fps[fpindex], "\n");
  return 0;
}

int diag_end(int fpindex) {
  return fclose(fps[fpindex]);
}

