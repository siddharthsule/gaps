#ifndef pdf_h_
#define pdf_h_

#include "LHAPDF/LHAPDF.h"

class pdf_wrapper {
  /**
   * @class pdf_wrapper
   * @brief Wrapper class for the PDF evaluation
   *
   * This class is a wrapper for the LHAPDF::PDF class. This provides
   * consistency with the GPU version, where a wrapper is needed for batch
   * evaluation. For CPU, this is a simple pass-through to LHAPDF methods.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables
  LHAPDF::PDF* lhapdf;

 public:
  // ---------------------------------------------------------------------------
  // constructor and destructor

  pdf_wrapper(const std::string& name = "CT14lo", int member = 0) {
    lhapdf = LHAPDF::mkPDF(name, member);
    LHAPDF::setVerbosity(0);
  }

  ~pdf_wrapper() { delete lhapdf; }

  // ---------------------------------------------------------------------------
  // member functions

  // Evaluate PDF for a single flavour, x and Q^2
  double xfxQ2(int fl, double x, double q2) const {
    return lhapdf->xfxQ2(fl, x, q2);
  }

  // Get the underlying LHAPDF::PDF pointer (for compatibility)
  LHAPDF::PDF* get_pdf() const { return lhapdf; }
};

#endif  // pdf_h_