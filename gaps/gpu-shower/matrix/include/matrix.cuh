#ifndef matrix_cuh_
#define matrix_cuh_

#include "event.cuh"
#include "interface.cuh"
#include "pdf.cuh"
#include "qcd.cuh"

class matrix {
  /**
   * @class matrix
   * @brief matrix element generation for e+ e- -> q qbar and p p -> Z
   */

 private:
  // ---------------------------------------------------------------------------
  // constants

 public:
  // Standard to Hard Event ratio
  double ws;

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
  // member functions

  // Matrix Element for e+e- -> qqbar, used for all LO
  __device__ double me2_ee2Zy2qq(int fl, double s, double t) const;

  // Matrix Element for q qbar to Z, for LHC NLO
  __device__ double me2qqZ(int fl, double s) const;
};

// -----------------------------------------------------------------------------
// Wrapper declarations

// function for unique process

// LEP LO
void lep_lo(thrust::device_vector<event>& d_events, matrix* matrix, int blocks,
            int threads);

// LEP NLO
void lep_nlo(thrust::device_vector<event>& d_events, matrix* matrix,
             alpha_s* as, int blocks, int threads);

// LHC LO
void lhc_lo(thrust::device_vector<event>& dv_events, matrix* matrix,
            pdf_wrapper* pdf, int blocks, int threads);

// LHC NLO
void lhc_nlo(thrust::device_vector<event>& d_events, matrix* matrix,
             alpha_s* as, pdf_wrapper* pdf, int blocks, int threads);

// all tasks wrapped in a function
void run_matrix(thrust::device_vector<event>& d_events, const params& p,
                int blocks);

#endif  // matrix_cuh_