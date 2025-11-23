#ifndef pdf_cuh_
#define pdf_cuh_

// Includes all modules, including LHAPDF
#include "base.cuh"

class pdf_wrapper {
  /**
   * @class pdf_wrapper
   * @brief Wrapper class for the PDF evaluation
   *
   * This class is a wrapper for the CuPDF class. This class is used to
   * evaluate the PDFs for a given set of flavours, x and Q^2 values.
   *
   * In the general event generation scenario, we cannot split events by flavour
   * - we need to evaluate the pdf for all events at once, for all flavours.
   * This means we evaluate pdf 11 times for each event, and then filter the
   * results to get the right value.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables
  CuPDF* pdf;

 public:
  // ---------------------------------------------------------------------------
  // constructor and destructor

  pdf_wrapper() : pdf(new CuPDF("CT14lo", 0)) { LHAPDF::setVerbosity(0); }
  ~pdf_wrapper() { delete pdf; }

  // ---------------------------------------------------------------------------
  // member functions

  // wrapper function that calls the kernel
  void evaluate(int* d_fl, double* d_x, double* d_q2, double* d_xf, int n,
                int blocks, int threads);
};

// After evaluating all flavours, filter data to retrieve right value
__global__ void select_flavour(int* fl, double* xf_in, double* xf_out, int n);

#endif  // pdf_cuh_