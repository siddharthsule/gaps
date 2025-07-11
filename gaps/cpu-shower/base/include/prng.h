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

// @brief generate a random number based on the seed
double generate_lcg(unsigned long& current_seed);

#endif  // prng_h_