#ifndef prng_cuh_
#define prng_cuh_

#include "base.cuh"

__device__ __inline__ double generate_lcg(unsigned long& current_seed) {
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

#endif  // prng_cuh_