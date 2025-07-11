#ifndef matrix_cuh
#define matrix_cuh

#include "event.cuh"
#include "qcd.cuh"

class matrix {
  /**
   * @class matrix
   * @brief matrix element generation
   *
   * this class is used to generate the leading order matrix element for the
   * ee->qq process the me^2 is calculated simulataneously
   * for all events, but with a few random numbers for flavour and direction.
   * this is a massless shower, so the system generates theoretical identical
   * events for all flavours.
   */

 private:
  // ---------------------------------------------------------------------------
  // constants

  double mz2, gz2, alpha, sin2tw;

 public:
  double amin, ye, ze, ws;

  //  LO/NLO and Energy
  bool nlo = false;
  double root_s = 0.;

  // Useful for Calculations
  double s = 0.;
  double s_hat = 0.;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  // constructor for device code
  __device__ void setup(bool nlo = false, double root_s = 91.2);

  // ---------------------------------------------------------------------------
  // get functions

  // ---------------------------------------------------------------------------
  // member functions

  // leading order matrix element generation
  __device__ double me2(int fl, double s, double t) const;
};

// -----------------------------------------------------------------------------
// kernels and wrappers

// cuda kernels to setup matrix and make lo points
// tip: cuda kernels cannot be member functions
__global__ void matrix_setup_kernel(matrix *matrix, bool nlo, double root_s);

__global__ void lep_lo(matrix *matrix, event *events, int n);
__global__ void lep_nlo(matrix *matrix, alpha_s *as, event *events, int n);

// all tasks wrapped in a function
void calc_lome(thrust::device_vector<event> &d_events, bool nlo, double root_s,
               double asmz, int blocks, int threads);

#endif  // matrix_cuh