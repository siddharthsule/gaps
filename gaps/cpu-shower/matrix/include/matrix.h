#ifndef matrix_h_
#define matrix_h_

#include "event.h"
#include "pdf.h"
#include "qcd.h"

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

  // PDF: defaults to NNPDF40MC (lo or nlo), or user-specified
  pdf_wrapper pdf;

  // Alpha S Calculator for CS NLO
  alpha_s as;

 public:
  // ---------------------------------------------------------------------------
  // constructor
  matrix(int process = 1, int nlo = 0, double root_s = 91.2,
         double asmz = 0.118, const std::string& pdf_name = "Null")
      : process(process),
        nlo(nlo),
        root_s(root_s),
        ws(0.25),
        pdf(pdf_name != "Null" ? pdf_name
            : nlo              ? "NNPDF40MC_nlo_as_01180"
                               : "NNPDF40MC_lo_as_01180"),
        as(asmz) {
    /**
     * @brief construct the matrix element generator
     *
     * stores the run configuration used by all matrix element functions: which
     * hard process to generate (lep or lhc), whether to run at lo or nlo, and
     * the collider center-of-mass energy. also builds the pdf and alpha_s(mz)
     * objects needed for the lhc and nlo calculations respectively.
     *
     * @param process the hard process: 1 = lep (e+e- -> qqbar), 2 = lhc (pp ->
     * Z)
     * @param nlo whether to generate at nlo (true) or lo (false)
     * @param root_s the collider center-of-mass energy
     * @param asmz the strong coupling alpha_s(mz) used to set up the alpha_s
     * runner
     * @param pdf_name the pdf set to use, or "Null" for the default NNPDF40MC
     * set
     */
  }

  // ---------------------------------------------------------------------------
  // member functions

  // Matrix Element for e+e- -> qqbar, used for all LO
  double me2_ee2Zy2qq(int fl, double s, double t) const;

  // Matrix Element for q qbar to Z, for LHC NLO
  double me2qqZ(int fl, double s) const;

  // function for unique process
  void lep_lo(event& ev);
  void lhc_lo(event& ev);
  void lep_nlo(event& ev);
  void lhc_nlo(event& ev);

  // wrapper for the matrix element
  void run(event& ev) {
    /**
     * @brief wrapper for matrix element calculation
     *
     * @param ev the event object
     */

    // LEP LO
    if ((process == 1) && !nlo) {
      lep_lo(ev);
    }

    // LEP NLO
    if ((process == 1) && nlo) {
      lep_nlo(ev);
    }

    // LHC LO
    if ((process == 2) && !nlo) {
      lhc_lo(ev);
    }

    // LHC NLO
    if ((process == 2) && nlo) {
      lhc_nlo(ev);
    }
  }
};

#endif  // matrix_h_