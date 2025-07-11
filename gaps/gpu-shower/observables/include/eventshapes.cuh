#ifndef eventshapes_cuh_
#define eventshapes_cuh_

#include "event.cuh"

// event shapes
__device__ void bubble_sort(vec4* moms, int n);
__global__ void calculate_ev_shapes(const event* events, double* results,
                                    int n);

#endif  // eventshapes_cuh_