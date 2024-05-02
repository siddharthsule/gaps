#ifndef jetrates_cuh_
#define jetrates_cuh_

#include "event.cuh"

// jet rates using durham algorithm
__device__ double yij(const vec4& p, const vec4& q, double ecm2);
__global__ void do_cluster(event* events, int n);

#endif  // jetrates_cuh_