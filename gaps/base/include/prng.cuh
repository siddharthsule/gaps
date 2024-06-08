#ifndef prng_cuh_
#define prng_cuh_

#include "base.cuh"

/**
 * PRNG with Linear congruential generator
 * ---------------------------------------
 *
 * To test that the two implementations are equivalent, we can use the random
 * number generator as the control variable. If the physics is correct, then
 * everthing should be the same.
 */

__device__ __inline__ double generate_lcg(unsigned long& current_seed) {
  // Update the seed
  current_seed = (lcg_a * (current_seed) + lcg_c) % lcg_m;

  // Return the random number
  return static_cast<double>(current_seed) / static_cast<double>(lcg_m);
}

#endif  // prng_cuh_