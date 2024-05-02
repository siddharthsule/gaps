#ifndef dalitz_cuh_
#define dalitz_cuh_

#include "event.cuh"

// dalitz plot
__global__ void calculate_dalitz(event* events, int n);

#endif  // dalitz_cuh_