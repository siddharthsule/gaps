#ifndef matrix_cuh
#define matrix_cuh

#include "event.cuh"
#include "pdf.cuh"
#include "qcd.cuh"

class matrix {
  /**
   * @class matrix
   * @brief matrix element generation
   *
   * this class is used to generate the leading order matrix element for the
   * ee->qq, eq->eq and qq-> ee processes the me^2 is calculated simulataneously
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

  // process, LO/NLO and Energy
  int process = 0;
  bool nlo = false;
  double root_s = 0.;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  // constructor for device code
  __device__ void setup(int process = 1, bool nlo = false,
                        double root_s = 91.2);

  // ---------------------------------------------------------------------------
  // get functions

  // ---------------------------------------------------------------------------
  // member functions

  // leading order matrix element generation
  __device__ double me2(int fl, double s, double t) const;

  // Matrix Element for q qbar to Z, for LHC NLO
  __device__ double me2qqZ(int fl, double s) const;
};

// -----------------------------------------------------------------------------
// kernels and wrappers

// cuda kernels to setup matrix and make lo points
// tip: cuda kernels cannot be member functions
__global__ void matrix_setup_kernel(matrix* matrix, int process, bool nlo,
                                    double root_s);

// LEP
__global__ void lep_lo(matrix* matrix, event* events, int n);
__global__ void lep_nlo(matrix* matrix, alpha_s* as, event* events, int n);

// LHC LO
void lhc_lo(thrust::device_vector<event>& dv_events, matrix* matrix, int blocks,
            int threads);

// LHC NLO
void lhc_nlo(thrust::device_vector<event>& d_events, matrix* matrix,
             alpha_s* as, int blocks, int threads);

// all tasks wrapped in a function
void calc_lome(thrust::device_vector<event>& d_events, int process, bool nlo,
               double root_s, double asmz, int blocks, int threads);

#endif  // matrix_cuh