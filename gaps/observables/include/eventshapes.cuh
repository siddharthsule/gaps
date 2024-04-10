#ifndef EVENTSHAPES_CUH_
#define EVENTSHAPES_CUH_

#include "event.cuh"

// Event Shapes
__device__ void bubbleSort(Vec4* moms, int n);
__global__ void calculateThr(Event* events, int N);
__global__ void calculateJetMBr(Event* events, int N);

#endif // EVENTSHAPES_CUH_