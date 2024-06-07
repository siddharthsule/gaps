#include "prng.h"

double generate_lcg(unsigned long& current_seed) {
  // Update the seed
  current_seed = (lcg_a * (current_seed) + lcg_c) % lcg_m;

  // Return the random number
  return static_cast<double>(current_seed) / static_cast<double>(lcg_m);
}