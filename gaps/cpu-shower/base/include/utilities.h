#ifndef utilities_h_
#define utilities_h_

#include "base.h"

// -----------------------------------------------------------------------------
// Small helpers with no home of their own

inline int choose_with_weights(const double* weights, int n, double rho) {
  /**
   * @brief Draw an index with probability proportional to its weight
   *
   * @param weights the weights, none of them negative
   * @param n how many there are
   * @param rho a random number in [0, 1)
   * @return the index drawn, or the last one where every weight is zero
   */

  double total = 0.;
  for (int i = 0; i < n; i++) total += weights[i];

  double target = rho * total;
  double running = 0.;
  for (int i = 0; i < n; i++) {
    running += weights[i];
    if (running > target) return i;
  }

  return n - 1;
}

#endif  // utilities_h_
