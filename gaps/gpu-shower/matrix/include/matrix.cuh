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
  // get functions

  // ---------------------------------------------------------------------------
  // member functions

  // leading order matrix element generation
  __device__ double me2_ee2Zy2qq(int fl, double s, double t) const;

  // Matrix Element for q qbar to Z, for LHC NLO
  __device__ double me2qqZ(int fl, double s) const;
};

// -----------------------------------------------------------------------------
// Wrapper declarations

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
void calc_lome(thrust::device_vector<event>& d_events, int process, bool nlo,
               double root_s, double asmz, int blocks, int threads,
               const std::string& pdf_name = "CT14lo");

#endif  // matrix_cuh