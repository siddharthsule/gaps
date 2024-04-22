#ifndef DALITZ_CUH_
#define DALITZ_CUH_

#include "event.cuh"

// Dalitz Plot
__global__ void calculateDalitz(Event* events, int N);

#endif  // DALITZ_CUH_