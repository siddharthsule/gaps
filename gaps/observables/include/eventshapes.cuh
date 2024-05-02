#ifndef eventshapes_cuh_
#define eventshapes_cuh_

#include "event.cuh"

// event shapes
__device__ void bubble_sort(vec4* moms, int n);
__global__ void calculate_thr(event* events, int n);
__global__ void calculate_jet_m_br(event* events, int n);

#endif  // eventshapes_cuh_