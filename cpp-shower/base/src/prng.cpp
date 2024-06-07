#include "prng.h"

void update_rng(unsigned long& current_seed, double& rand) {
  current_seed = (lcg_a * (current_seed) + lcg_c) % lcg_m;
  rand = static_cast<double>(current_seed) / lcg_m;
}