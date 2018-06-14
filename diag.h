#pragma once

int diag_init(char* filename);
int diag_write(int fpindex, double x, double y);
int diag_end(int fpindex);
int diag_write_list(int fpindex, double x, double* y, int size);
