#ifndef jets_cuh_
#define jets_cuh_

#include "event.cuh"

// jet rates using durham algorithm
__device__ double yij(const vec4& p, const vec4& q, double ecm2);
__global__ void cluster_durham(const event* events, double* results, int n);

#endif  // jets_cuh_