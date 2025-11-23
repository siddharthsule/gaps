#ifndef qcd_cuh_
#define qcd_cuh_

#include "event.cuh"

// QCD constants
const double k_tr = 0.5;
const double k_ca = k_nc;
const double k_cf = (k_nc * k_nc - 1.) / (2. * k_nc);

class alpha_s {
  /**
   * @class alpha_s
   * @brief The strong coupling constant
   *
   * This class is used to calculate the strong coupling constant at a given
   * scale.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  int n_loops;
  double mc2, mb2, mz2, asmz, asmb, asmc;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  __device__ alpha_s(double mz, double asmz, int n_loops = 2, double mb = 4.75,
                     double mc = 1.3);

  // ---------------------------------------------------------------------------
  // member functions

  // setup function for device code
  __device__ void setup(double mz, double asmz, int n_loops = 2,
                        double mb = 4.75, double mc = 1.3);

  // all the required functions to calculate the strong coupling constant
  __device__ double beta0(int nf) const;
  __device__ double beta1(int nf) const;
  __device__ double as1(double t) const;
  __device__ double as2(double t) const;
  __device__ double operator()(double t) const;
};

// setup the alpha_s class
__global__ void as_setup_kernel(alpha_s* as, double mz, double asmz,
                                int n_loops = 2, double mb = 4.75,
                                double mc = 1.27);

// calculate the strong coupling constant
__global__ void as_value(alpha_s* as, double* as_val, double t);

#endif  // qcd_cuh_