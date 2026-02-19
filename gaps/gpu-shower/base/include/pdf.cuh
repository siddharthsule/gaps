#ifndef pdf_cuh_
#define pdf_cuh_

#include "base.cuh"
#include "event.cuh"
#include <string>

/**
 * How to avoid the header-only nature of LHAPDF from stopping us? We can
 * define the class here, and implement everything in the .cu file. This
 * way, the PDF evaluation kernel is written to device once only!
 *
*/
class CuPDF;
namespace LHAPDF {
  void setVerbosity(int);
}

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

  // Define the constructor in pdf.cu
  pdf_wrapper(const std::string& name = "CT14lo", int member = 0);
  ~pdf_wrapper();

  // ---------------------------------------------------------------------------
  // member functions

  // wrapper function that calls the kernel
  void evaluate(int* d_fl, double* d_x, double* d_q2, double* d_xf, int n,
                int blocks, int threads);

  // Multiply cross section by PDF values
  void multiply_dxs_by_pdf(event* events, int n, double* xf, int blocks,
                           int threads);
};




#endif  // pdf_cuh_