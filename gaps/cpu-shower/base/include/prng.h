#ifndef prng_h_
#define prng_h_

#include "base.h"

/**
 * PRNG with Linear congruential generator
 * ---------------------------------------
 *
 * To test that the two implementations are equivalent, we can use the random
 * number generator as the control variable. If the physics is correct, then
 * everthing should be the same.
 */

inline double generate_lcg(unsigned long& current_seed) {
  /**
   * @brief Generate a random number using a Linear Congruential Generator (LCG)
   *
   * @param current_seed The current seed for the LCG
   * @return The random number
   */

  // Update the seed
  current_seed = (lcg_a * (current_seed) + lcg_c) % lcg_m;

  // Return the random number
  return static_cast<double>(current_seed) / static_cast<double>(lcg_m);
}

#endif  // prng_h_