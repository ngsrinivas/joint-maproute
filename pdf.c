#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NBINS 10001
#define BUFSIZE 10000

/*
Requirements: shouldn't have to care about number of numbers
shouldn't have to care about range of numbers
It's okay to do two-pass of array but not of the file
Be able to input floating point numbers.
*/

struct numbuffer {
  double nums[BUFSIZE];
  struct numbuffer* next;
};

double pdf[NBINS];

int print_all_same(double num) {
  
  int i;
  for (i=0; i<=100; i++) {
    printf("%lf %lf\n", num, i/100.0);
  }
  
  return 0;
  
}

int main(int argc, char** argv) {

  // first pass: store numbers read, get the range of numbers entered
  struct numbuffer* first;
  first = (struct numbuffer*)malloc(sizeof(struct numbuffer));
  first->next = NULL;
  
  double thisnum = 0.0;
  int num_nums = 0; // Number of numbers filled in current buffer
  double min_num = 0.0, max_num = 0.0; // global min and max numbers, useful for binning later on
  struct numbuffer* curr_buffer = first;
  int total_nums = 0; // total number of numbers
  int buffer_count = 0;

  fprintf(stderr, "Reading in from input..\n");

  while(scanf("%lf\n", &thisnum) > 0) { // first scan of numbers
    num_nums ++;
    total_nums ++;
    
    if (num_nums > BUFSIZE) {
      buffer_count ++;
      fprintf(stderr, "Read in one buffer full.. buffer %d\n", buffer_count);
      curr_buffer->next = (struct numbuffer*)malloc(sizeof(struct numbuffer));
      curr_buffer = curr_buffer->next;
      curr_buffer->next = NULL; // set null to indicate end of the list
      num_nums = 1;
    }

    curr_buffer->nums[num_nums - 1] = thisnum;
    
    // adjust max and min
    if (total_nums > 1) {
      if (thisnum > max_num) max_num = thisnum;
      if (thisnum < min_num) min_num = thisnum;
    }
    else max_num = min_num = thisnum;
    
  }

  int last_index;

  double multfactor = (NBINS-1)/(max_num - min_num);

  // check if there is a range in the distribution
  // if not, just print this number for all percentiles.
  if (fabs(max_num - min_num) <= 0.000001) {
    fprintf(stderr, "Range is quite small.\n");
    return print_all_same(max_num);
  }
  
  int i, j;
  curr_buffer = first;
  
  fprintf(stderr, "Finished reading numbers: last buffer %d is %d full\n", buffer_count + 1, num_nums);
  fprintf(stderr, "Starting to bin\n");

  buffer_count = 0;

  while (curr_buffer != NULL) { // start binning items from buffers into bins
    if (curr_buffer->next == NULL) last_index = num_nums;
    else last_index = BUFSIZE;

    for (i=0; i<last_index; i++)
      pdf[(int)((curr_buffer->nums[i] - min_num) * multfactor)] += 1.0;
    
    curr_buffer = curr_buffer->next;
    buffer_count ++;
    fprintf(stderr, "Binned buffer %d, moving to next..\n", buffer_count);
  }

  fprintf(stderr, "Normalizing\n");

  // normalize bins
  for (i=0; i<NBINS; i++) {
    pdf[i] /= total_nums;
  }
  
  fprintf(stderr, "Printing CDF values\n");

  double cdf = 0.0; // for cumulative values
  
  // print pdf value with bin starting value
  for (i=0; i<NBINS; i++) {
    cdf += pdf[i];
    printf("%lf %lf\n", min_num + (i/multfactor), cdf);
  }

  return 0;
}

