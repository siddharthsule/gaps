#ifndef QCD_CUH_
#define QCD_CUH_

// Parton includes Base, which has the CUDA libraries
#include "parton.cuh"

/**
 * The Strong Coupling Constant
 * ----------------------------
 *
 * This file contains the QCD constants and the alpha_s class. The alpha_s class
 * is a simple class that calculates the strong coupling constant at a given
 * scale. The class is designed to be used in a CUDA kernel, so it is a simple
 * class with no dynamic memory allocation.
 */

// QCD Constants, maybe you can use for SU(Not 3) ?
const double kNC = 3.0;
const double kTR = 0.5;
const double kCA = kNC;
const double kCF = (kNC * kNC - 1.0) / (2.0 * kNC);

// A lot of the member functions can be defined here, but we need a .cu file to
// define the kernels, so might as well define everything in the .cu file!
class AlphaS {
 private:
  int order;
  double mc2, mb2, mz2, asmz, asmb, asmc;

 public:
  // Constructor
  __device__ AlphaS(double mz, double asmz, int order = 1, double mb = 4.75,
                    double mc = 1.27);

  // Setup function for device code
  __device__ void setup(double mz, double asmz, int order = 1, double mb = 4.75,
                        double mc = 1.27);

  // All the required functions to calculate the strong coupling constant
  __device__ double Beta0(int nf);
  __device__ double Beta1(int nf);
  __device__ double As0(double t);
  __device__ double As1(double t);
  __device__ double operator()(double t);
};

// Setup the alpha_s class
__global__ void asSetupKernel(AlphaS *as, double mz, double asmz, int order = 1,
                              double mb = 4.75, double mc = 1.27);

// Calculate the strong coupling constant
__global__ void asValue(AlphaS *as, double *asval, double t);
__global__ void asShowerKernel(AlphaS *as, Event *events, int N);

#endif  // QCD_CUH_