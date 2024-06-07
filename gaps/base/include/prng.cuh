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

__device__ __inline__ void update_rng(unsigned long& current_seed,
                                      double& rand) {
  current_seed = (lcg_a * (current_seed) + lcg_c) % lcg_m;
  rand = static_cast<double>(current_seed) / lcg_m;
}

#endif  // prng_cuh_