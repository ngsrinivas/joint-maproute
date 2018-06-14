#include<stdio.h>
#include<stdlib.h>
#include "params.h"

// speed of transmission through the backbone
// in miles per millisecond
#define LIGHT_SPEED 130.39767

// things to capture. Interfaces: through 
// multiple files whose file descriptors are provided.

// general inputs / "basic" inputs
// Note: The previous input file for starting off the "random instance" simulations
// CANNOT BE USED HERE for reading through this function -- the parameters used here
// are in a slightly different order and there are different parameters present 
// and absent here.
int read_basic_inputs(FILE* fp, struct params_t* p) {
  fscanf(fp, "%d\n", &p->I);
  fscanf(fp, "%d\n", &p->J);
  fscanf(fp, "%d\n", &p->N);
  fscanf(fp, "%d\n", &p->C);
  fscanf(fp, "%lf\n", &p->K);
  fscanf(fp, "%lf\n", &p->L); // coefficient for cost. Can be set to zero for performance only optimization.
  return 0;
}

// number of links on every DC
int read_link_information(FILE* fp, struct params_t* p) {
  int i;
  int total = 0;
  for (i=0; i<p->I; i++) {
    fscanf(fp, "%d\n", &p->lnum[i]);
    total += p->lnum[i];
  }
  p->total_links = total; // assign total
  return 0;
}

// traffic volumes by prefix (there's a time series)
int read_traffic_volumes(FILE* fp, struct params_t *p, double* timestamp) {
  fscanf(fp, "%lf ", timestamp);
  int c, n;
  double TT = 0.0;
  for (c=0; c<p->C; c++) {
    for (n=0; n<p->N; n++) {
      fscanf(fp, "%lf ", &p->vol[c*p->N + n]);
      TT += p->vol[c*p->N + n];
    }
  }
  p->TT = TT; // assign total traffic amount
  fscanf(fp, "\n");
  return 0;
}

// performance by prefix (there's time series)
int read_performance(FILE* fp, struct params_t* p, double* timestamp) {
  *timestamp = 0.0;
  int I, J, i, j, c, k, j_i;
  int C = p->C;
  
  fscanf(fp, "%d %d\n", &I, &J);
  
  for (i=0; i<I; i++) {
    fscanf(fp, "%d\n", &j_i);
    for (j=0; j<j_i; j++) { // each line of the input
      for (c=0; c<C; c++) {    
	fscanf(fp, "%lf ", &p->perf[i*J*C + j*C + c]);
      }
      fscanf(fp, "\n");
    }
  }
  
  return 0;
  
}

// link capacities
int read_link_capacities(FILE* fp, struct params_t* p) {
  int i, j, k, I, J, j_i;
  
  fscanf(fp, "%d %d\n", &I, &J);
  for (i=0; i<I; i++) { // each line of input
    fscanf(fp, "%d ", &j_i);
    for (j=0; j<j_i; j++) {
      fscanf(fp, "%lf ", &p->cap[i*J + j]);
    }
    fscanf(fp, "\n");
  }
  return 0;
}

// link pricing
int read_link_pricing(FILE* fp, struct params_t* p) {
  int i, j, k, I, J, j_i;
  
  fscanf(fp, "%d %d\n", &I, &J);
  
  for (i=0; i<I; i++) { // for each line of input
    fscanf(fp, "%d ", &j_i);
    for (j=0; j<j_i; j++) {
      fscanf(fp, "%lf ", &p->price[i*J + j]);
    }
    fscanf(fp, "\n");
  }
  return 0;
}

int read_w_eps(FILE* fp, struct params_t* p) {
  // This function does not read load balancing parameters from anywhere.
  // Since we are not evaluating this in anyway, for now I will just set
  // them so that those constraints are always feasible, and hence do not 
  // affect the problem.
  int i;
  for (i=0; i<p->I; i++) {
    p->w[i] = 0.5;
    p->eps[i] = 1.0;
  }
  return 0;
}
    
// Input master for one instance. Read one instance of the
// problem. No assumptions on time series.
// "files" is the list of string filenames for one time access to read
// a single instance.
// This function takes care of allocating an appropriate amount of
// memory for parameters it reads in, and returns a filled-in
// parameter structure.

// IMPORTANT NOTE: For different experiments, use differnet
// input_master functions. The assumptions these functions make on the
// data received and the way they handle the file descriptors are
// (necessarily) very different.

int input_master_single_instance(char** files, struct params_t* p) {
  
  // ordering on files:
  // 1. basic inputs
  // 2. link number information
  // 3. traffic information
  // 4. performance information
  // 5. capacity information
  // 6. pricing information

  FILE* fbasic_inputs = fopen(files[0], "r");
  FILE* flinknum_info = fopen(files[1], "r");
  FILE* ftraffic = fopen(files[2], "r");
  FILE* fperf = fopen(files[3], "r");
  FILE* fcaps = fopen(files[4], "r");
  FILE* fpricing = fopen(files[5], "r");

  if (fbasic_inputs == NULL || flinknum_info == NULL || ftraffic == NULL || fperf == NULL || fcaps == NULL || fpricing == NULL) {
    fprintf(stderr, "Inputs: Couldn't open files for reading inputs for simulation!\n");
    fprintf(stderr, "%s %s %s %s %s %s\n", files[0], files[1], files[2], files[3], files[4], files[5]);
    exit(0);
  }

  // for now, assume that the FPs can be directly fed to the
  // input reading functions, and setup the problem parameters
  // structure p appropriately.

  // Read in basic inputs which dictate memory allocation for the rest.
  read_basic_inputs(fbasic_inputs, p);
  
  // setup link count information for DCs
  p->lnum = (int*)malloc(p->I * sizeof(int));
  if (p->lnum == NULL) { fprintf(stderr, "Inputs: Memory allocation for link counts failed!\n"); exit(0); }
  read_link_information(flinknum_info, p);
  
  double timestamp;
  
  // setup traffic information
  p->vol = (double*)malloc(p->C * p->N * sizeof(double));
  if (p->vol == NULL) { fprintf(stderr, "Inputs: Memory allocation for traffic volumes failed!\n"); exit(0); }
  int ninstances;
  fscanf(ftraffic, "%d\n", &ninstances);
  read_traffic_volumes(ftraffic, p, &timestamp);
  
  // setup performance information
  p->perf = (double*)malloc(p->I * p->J * p->C * sizeof(double));
  if (p->perf == NULL) { fprintf(stderr, "Inputs: Memory allocation for performance failed!\n"); exit(0); }
  read_performance(fperf, p, &timestamp);

  // setup capacity information
  p->cap = (double*)malloc(p->I * p->J * sizeof(double));
  if (p->cap == NULL) { fprintf(stderr, "Inputs: Memory allocation for capacity failed!\n"); exit(0); }
  read_link_capacities(fcaps, p);
  
  // setup pricing information
  p->price = (double*)malloc(p->I * p->J * sizeof(double));
  if (p->price == NULL) { fprintf(stderr, "Inputs: Memory allocation for link pricing failed!\n"); exit(0); }
  read_link_pricing(fpricing, p);
  
  p->w = (double*)malloc(p->I * sizeof(double));
  p->eps = (double*)malloc(p->I * sizeof(double));
  if (p->w == NULL || p->eps == NULL) { fprintf(stderr, "Inputs: Memory allocation for DC load balancing failed!\n"); exit(0); }
  read_w_eps(NULL, p);

  // close all open file descriptors
  fclose(fbasic_inputs);
  fclose(flinknum_info);
  fclose(ftraffic);
  fclose(fperf);
  fclose(fcaps);
  fclose(fpricing);

  return 0;
  
}

// Deallocation master for single instance problem inputting.
// frees all memory from a single params_t structure.
int dealloc_master_single_instance(struct params_t* p) {
  free(p->price);
  free(p->perf);
  free(p->vol);
  free(p->w);
  free(p->eps);
  free(p->cap);
  free(p->lnum);
  return 0;
}

// Input master for traffic series. Input parameters include the same
// set as the single instance master, but here every time the function
// is called you get a new instance from the traffic series. An
// additional parameter informs you of number of instances involved.
// Just call read_traffic_volumes() ninstances-1 times to get
// the rest of the traffic series in the params_t object. For this
// purpose, the traffic file descriptor is also saved and returned in
// one of the arguments.

int input_master_traffic_timeseries(char** files, struct params_t* p, int* ninstances, FILE** ftraffic_save) {
  
  // ordering on files:
  // 1. basic inputs
  // 2. link number information
  // 3. traffic information
  // 4. performance information
  // 5. capacity information
  // 6. pricing information

  FILE* fbasic_inputs = fopen(files[0], "r");
  FILE* flinknum_info = fopen(files[1], "r");
  FILE* ftraffic = fopen(files[2], "r");
  FILE* fperf = fopen(files[3], "r");
  FILE* fcaps = fopen(files[4], "r");
  FILE* fpricing = fopen(files[5], "r");

  if (fbasic_inputs == NULL || flinknum_info == NULL || ftraffic == NULL || fperf == NULL || fcaps == NULL || fpricing == NULL) {
    fprintf(stderr, "Inputs: Couldn't open files for reading inputs for simulation!\n");
    fprintf(stderr, "%s %s %s %s %s %s\n", files[0], files[1], files[2], files[3], files[4], files[5]);
    exit(0);
  }
  
  // Save ftraffic for future traffic series reads.
  *ftraffic_save = ftraffic;

  // for now, assume that the FPs can be directly fed to the
  // input reading functions, and setup the problem parameters
  // structure p appropriately.

  // Read in basic inputs which dictate memory allocation for the rest.
  read_basic_inputs(fbasic_inputs, p);
  
  // setup link count information for DCs
  p->lnum = (int*)malloc(p->I * sizeof(int));
  if (p->lnum == NULL) { fprintf(stderr, "Inputs: Memory allocation for link counts failed!\n"); exit(0); }
  read_link_information(flinknum_info, p);
  
  double timestamp;
  
  // setup traffic information
  p->vol = (double*)malloc(p->C * p->N * sizeof(double));
  if (p->vol == NULL) { fprintf(stderr, "Inputs: Memory allocation for traffic volumes failed!\n"); exit(0); }
  fscanf(ftraffic, "%d\n", ninstances);
  read_traffic_volumes(ftraffic, p, &timestamp);
  
  // setup performance information
  p->perf = (double*)malloc(p->I * p->J * p->C * sizeof(double));
  if (p->perf == NULL) { fprintf(stderr, "Inputs: Memory allocation for performance failed!\n"); exit(0); }
  read_performance(fperf, p, &timestamp);

  // setup capacity information
  p->cap = (double*)malloc(p->I * p->J * sizeof(double));
  if (p->cap == NULL) { fprintf(stderr, "Inputs: Memory allocation for capacity failed!\n"); exit(0); }
  read_link_capacities(fcaps, p);
  
  // setup pricing information
  p->price = (double*)malloc(p->I * p->J * sizeof(double));
  if (p->price == NULL) { fprintf(stderr, "Inputs: Memory allocation for link pricing failed!\n"); exit(0); }
  read_link_pricing(fpricing, p);
  
  p->w = (double*)malloc(p->I * sizeof(double));
  p->eps = (double*)malloc(p->I * sizeof(double));
  if (p->w == NULL || p->eps == NULL) { fprintf(stderr, "Inputs: Memory allocation for DC load balancing failed!\n"); exit(0); }
  read_w_eps(NULL, p);

  // close all open file descriptors
  fclose(fbasic_inputs);
  fclose(flinknum_info);
  // Can't close ftraffic because it'll be used later!
  // fclose(ftraffic);
  fclose(fperf);
  fclose(fcaps);
  fclose(fpricing);

  return 0;
  
}

// Read distance-matrix file and print into provided backbone 
// latency array
int read_backbone_lats_from_distances(char* distfile, double* blats) {
  
  FILE* fp = fopen(distfile, "r");
  if (fp == NULL) {
    printf("Could not open backbone link distances file! %s\n", distfile);
    exit(0);
  }
  
  int I;
  int i, e;
  fscanf(fp, "%d\n", &I);
  for (i=0; i<I-1; i++) {
    for (e=i+1; e<I; e++) {
      fscanf(fp, "%lf ", &blats[i*I + e]);
      blats[i*I + e] /= LIGHT_SPEED;
      blats[e*I + i] = blats[i*I + e];
    }
    fscanf(fp, "\n");
  }

  fclose(fp);
  
  return 0;
  
}
