#ifndef matrix_h_
#define matrix_h_

#include "event.h"
#include "qcd.h"

class matrix {
  /**
   * @class matrix
   * @brief matrix element generation
   *
   * this class is used to generate the leading order matrix element for the
   * ee->qq process the me^2 is calculated simulataneously for all events,
   * but with a few random numbers for flavour and direction. this is a
   * massless shower, so the system generates theoretical identical
   * events for all flavours.
   */

 private:
  // ---------------------------------------------------------------------------
  // constants

  double mz2, gz2, alpha, sin2tw;

 public:
  double amin, ye, ze, ws;

  // LO/NLO and Energy
  bool nlo = false;
  double root_s = 0.;

  // Useful for Calculations
  double s = 0.;
  double s_hat = 0.;

  // Alpha S Calculator for CS NLO
  alpha_s as;

 public:
  // ---------------------------------------------------------------------------
  // constructor
  matrix(bool nlo = false, double root_s = 91.2, double asmz = 0.118)
      : nlo(nlo),
        root_s(root_s),
        s(pow(root_s, 2.)),
        s_hat(s),
        mz2(pow(mz, 2.)),
        gz2(pow(gz, 2.)),
        alpha(1. / 128.802),
        sin2tw(0.22293),
        amin(1.e-10),
        ws(0.25),
        as(mz, asmz) {}

  // ---------------------------------------------------------------------------
  // member functions

  // Matrix Element for e+e- -> qqbar, used for all LO
  double me2(int fl, double s, double t) const;

  // function for unique processes
  void lep_lo(event& ev);
  void lep_nlo(event& ev);

  // wrapper for the matrix element
  void run(event& ev) {
    // LEP LO
    if (!nlo) {
      lep_lo(ev);
    }

    // LEP NLO
    else {
      lep_nlo(ev);
    }
  }
};

#endif  // matrix_h_