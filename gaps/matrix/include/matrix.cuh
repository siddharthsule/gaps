#ifndef matrix_cuh
#define matrix_cuh

// parton includes base, which has the cuda libraries
#include "event.cuh"

/**
 * matrix element generation
 * -------------------------
 *
 * this class is used to generate the leading order matrix element for the
 * e+e- -> q qbar process. the me^2 is calculated simulataneously for all
 * events, but with a few random numbers for flavour and direction. this is a
 * massless shower, so the system generates theoretical identical events for
 * all flavours.
 */
class matrix {
 private:
  double alphas, ecms, mz2, gz2, alpha, sin2tw, amin, ye, ze, ws;

 public:
  // constructor
  matrix(double alphas = asmz, double ecms = 91.2);

  // setup function for device code
  __device__ void setup(double alphas = asmz, double ecms = 91.2);

  // leading order matrix element generation
  __device__ double me2(int fl, double s, double t);

  // getters
  __device__ double get_ecms() { return ecms; };
};

// cuda kernels to setup matrix and make lo points
// tip: cuda kernels cannot be member functions
__global__ void matrix_setup_kernel(matrix* matrix, double e);
__global__ void lo_point_kernel(matrix* matrix, event* ev, int n);

// all tasks wrapped in a function
void calc_lome(thrust::device_vector<event>& d_events, double e);

#endif  // matrix_cuh