#ifndef JETRATES_CUH_
#define JETRATES_CUH_

#include "parton.cuh"

// Jet rates using Durham algorithm
__device__ double Yij(const Vec4& p, const Vec4& q, double ecm2);
__global__ void doCluster(Event* events, int N);


#endif // JETRATES_CUH_