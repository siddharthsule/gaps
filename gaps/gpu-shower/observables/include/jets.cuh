#ifndef jets_cuh_
#define jets_cuh_

#include "event.cuh"

// jet rates using durham algorithm
__device__ double yij(const vec4& p, const vec4& q, double ecm2);
__global__ void cluster_durham(const event* events, double* results, int n);

// jet clustering using the generalized kt algorithm
__device__ void bubble_sort_pt(vec4* moms, int n);
__device__ double dR2(const vec4& p, const vec4& q);
__device__ double dij(const vec4& p, const vec4& q);
__global__ void cluster_genkt(const event* events, double* results, int n);

#endif  // jets_cuh_