#ifndef qcd_cuh_
#define qcd_cuh_

// parton includes base, which has the cuda libraries
#include "event.cuh"

/**
 * the strong coupling constant
 * ----------------------------
 *
 * this file contains the qcd constants and the alpha_s class. the alpha_s class
 * is a simple class that calculates the strong coupling constant at a given
 * scale. the class is designed to be used in a cuda kernel, so it is a simple
 * class with no dynamic memory allocation.
 */

// qcd constants, maybe you can use for su(not 3) ?
const double k_nc = 3.;
const double k_tr = 0.5;
const double k_ca = k_nc;
const double k_cf = (k_nc * k_nc - 1.) / (2. * k_nc);

// a lot of the member functions can be defined here, but we need a .cu file to
// define the kernels, so might as well define everything in the .cu file!
class alpha_s {
 private:
  int order;
  double mc2, mb2, mz2, asmz, asmb, asmc;

 public:
  // constructor
  __device__ alpha_s(double mz, double asmz, int order = 1, double mb = 4.75,
                     double mc = 1.27);

  // setup function for device code
  __device__ void setup(double mz, double asmz, int order = 1, double mb = 4.75,
                        double mc = 1.27);

  // all the required functions to calculate the strong coupling constant
  __device__ double beta0(int nf);
  __device__ double beta1(int nf);
  __device__ double as0(double t);
  __device__ double as1(double t);
  __device__ double operator()(double t);
};

// setup the alpha_s class
__global__ void as_setup_kernel(alpha_s *as, double mz, double asmz,
                                int order = 1, double mb = 4.75,
                                double mc = 1.27);

// calculate the strong coupling constant
__global__ void as_value(alpha_s *as, double *asval, double t);
__global__ void as_shower_kernel(alpha_s *as, event *events, double *asval,
                                 int n);

#endif  // qcd_cuh_