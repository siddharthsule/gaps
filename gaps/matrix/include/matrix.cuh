#ifndef MATRIX_CUH
#define MATRIX_CUH

// Parton includes Base, which has the CUDA libraries
#include "event.cuh"

/**
 * Matrix Element Generation
 * -------------------------
 *
 * This class is used to generate the leading order matrix element for the
 * e+e- -> q qbar process. The ME^2 is calculated simulataneously for all
 * events, but with a few random numbers for flavour and direction. This is a
 * massless shower, so the system generates theoretical identical events for
 * all flavours.
 */
class Matrix {
 private:
  double alphas, ecms, MZ2, GZ2, alpha, sin2tw, amin, ye, ze, ws;

 public:
  // Constructor
  Matrix(double alphas = asmz, double ecms = 91.2);

  // Setup function for device code
  __device__ void setup(double alphas = asmz, double ecms = 91.2);

  // Leading Order Matrix Element Generation
  __device__ double ME2(int fl, double s, double t);

  // Getters
  __device__ double GetECMS() { return ecms; };
};

// CUDA Kernels to Setup Matrix and make LO Points
// TIP: CUDA KERNELS CANNOT BE MEMBER FUNCTIONS
__global__ void matrixSetupKernel(Matrix* matrix, double E);
__global__ void loPointKernel(Matrix* matrix, Event* ev, int N);

// All tasks wrapped in a function
void calcLOME(thrust::device_vector<Event>& d_events, double E);

#endif  // MATRIX_CUH