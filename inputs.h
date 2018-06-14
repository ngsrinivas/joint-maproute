#pragma once

#include<stdlib.h>

#include "params.h"

int input_master_single_instance(char**, struct params_t*);
int input_master_traffic_timeseries(char**, struct params_t*, int*, FILE**);
int read_backbone_lats_from_distances(char*, double*);
int read_traffic_volumes(FILE*, struct params_t*, double*);
