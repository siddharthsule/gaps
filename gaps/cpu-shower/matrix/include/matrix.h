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

  // PDF: CT14lo
  pdf_wrapper pdf;

  // Alpha S Calculator for CS NLO
  alpha_s as;

 public:
  // ---------------------------------------------------------------------------
  // constructor
  matrix(int process = 1, bool nlo = false, double root_s = 91.2,
         double asmz = 0.118, const std::string& pdf_name = "CT14lo")
      : process(process),
        nlo(nlo),
        root_s(root_s),
        ws(0.25),
        pdf(pdf_name),
        as(asmz) {}

  // ---------------------------------------------------------------------------
  // member functions

  // Matrix Element for e+e- -> qqbar, used for all LO
  double me2_ee2Zy2qq(int fl, double s, double t) const;

  // Matrix Element for q qbar to Z, for LHC NLO
  // double me2qqgamma(int fl, double s);
  double me2qqZ(int fl, double s) const;

  // function for unique process
  void lep_lo(event& ev);
  void lhc_lo(event& ev);
  void lep_nlo(event& ev);
  void lhc_nlo(event& ev);

  // wrapper for the matrix element
  void run(event& ev) {
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