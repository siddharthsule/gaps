#ifndef MATRIX_CUH
#define MATRIX_CUH

#include <cmath>

#include "particle.cuh"

class Matrix {
public:

  // Constructor
  Matrix(double alphas, double ecms = 91.2);

  // Setup function for device code
  __device__ void setup(double alphas, double ecms);

  // Leading Order Matrix Element Generation
  __device__ double ME2(int fl, double s, double t);

private:
  double alphas, ecms, MZ2, GZ2, alpha, sin2tw, amin, ye, ze, ws;
};

// Declaration of the CUDA kernels, CANNOT BE MEMBER FUNCTIONS
__global__ void matrixSetupKernel(Matrix* matrix, double alphas, double ecms);
__global__ void loPointKernel(Matrix* matrix, Event* ev, int n);

#endif // MATRIX_CUH
