#include "mosek.h"

/*** Helper functions ****/

// Some useful definitions.
// link load coefficients for constraints
int y_coeff[6] = {1, 3, 10, 70, 5000, 5000};

// Capacity coefficients
double cap_coeff[6] = {0.0, 2.0/3.0, 16.0/3.0, 178.0/3.0, 14968.0/3.0, 14968.0/3.0};



// compute total amount of traffic, given vol
// per client per mapping node
double total_traffic(double* vol, int size) {
  int i;
  double tt = 0.0;
  for (i=0; i<size; i++) {
    tt += vol[i];
  }
  return tt;
}

/* int main() { */
/*   printf("Hello, world!\n"); */
/*   void *x; */
/*   printstr(x, "foo\n"); */
/*   return 0; */
/* } */

// prints MOSEK output to the terminal
static void MSKAPI printstr(void *handle, char str[]) {
  printf("%s", str);
}

// initialize a mosek environment
int init_mosek_env(MSKenv_t* env) {
  MSKrescodee r;
  
  r = MSK_makeenv(env, NULL, NULL, NULL, NULL);
  
  /* if (r == MSK_RES_OK) */
  /*   MSK_linkfunctoenvstream(env, MSK_STREAM_LOG, NULL, printstr); */

  if (r == MSK_RES_OK)
    r = MSK_initenv(*env);
  
  if (r == MSK_RES_OK) return 0;
  else {
    printf("Could not initialize environment properly!\n");
    exit(0);
  }

}
