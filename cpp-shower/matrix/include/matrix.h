#ifndef matrix_h_
#define matrix_h_

// parton includes base, which has the cuda libraries
#include "event.h"

/**
 * matrix element generation
 * -------------------------
 *
 * this class is used to generate the leading order matrix element for the
 * e+e- -> q qbar process. the me^2 is calculated simulataneously for all
 * events, but with a few random numbers for flavour and direction. this is a
 * massless shower, so the system generates theoretical identical events for
 * all flavours.
 */

class matrix {
 private:
  double alphas, ecms, mz2, gz2, alpha, sin2tw, amin, ye, ze, ws;

 public:
  // constructor
  matrix(double alphas = asmz, double ecms = 91.2);

  // leading order matrix element generation
  double me2(int fl, double s, double t);

  // generate a leading order point
  void generate_lo_point(event &ev);
};

#endif  // matrix_h_