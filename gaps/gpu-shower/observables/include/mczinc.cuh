#ifndef mczinc_cuh_
#define mczinc_cuh_

#include "event.cuh"

// mc zinc variables
__global__ void calculate_mczinc(const event* events, double* results, int n);

#endif  // mczinc_cuh_