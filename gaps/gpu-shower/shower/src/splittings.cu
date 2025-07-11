#include "shower.cuh"

/**
 * Here you can define and edit the splitting functions for the shower. Make
 * sure to update the estimates, integrals and generate_z functions as well,
 * or else the vetoalgorithm will not work correctly.
 *
 * the functions are defined as (type)_(splitting)_(case):
 * - Type: V = Value, E = Estimate, I = Integral, G = Generate Z
 * - Splitting: q -> qg, q -> gq, g -> gg, g -> qq
 * - Case: FF, FI, IF, II (Not in GAPS 1)
 *
 * To understand how these are Implemented, see the code after the functions.
 */

// -----------------------------------------------------------------------------
// FF splittings: z and y

// Values --------------------------------

__device__ double v_qqg_ff(double z, double y) {
  return k_cf * (2. / (1. - z * (1. - y)) - (1. + z));
}

__device__ double v_ggg_ff(double z, double y) {
  return k_ca / 2. * (2. / (1. - z * (1. - y)) - 2. + z * (1. - z));
}

__device__ double v_gqq_ff(double z, double y) {
  return k_tr / 2. * (1. - 2. * z * (1. - z));
}

// Estimates ------------------------------

__device__ double e_qqg_ff(double z) { return k_cf * 2. / (1. - z); }

__device__ double e_ggg_ff(double z) { return k_ca / (1. - z); }

__device__ double e_gqq_ff(double z) { return k_tr / 2.; }

// Integrals ------------------------------

__device__ double i_qqg_ff(double zm, double zp) {
  return k_cf * 2. * log((1. - zm) / (1. - zp));
}

__device__ double i_ggg_ff(double zm, double zp) {
  return k_ca * log((1. - zm) / (1. - zp));
}

__device__ double i_gqq_ff(double zm, double zp) {
  return k_tr / 2. * (zp - zm);
}

// Generate Z ------------------------------

__device__ double g_qqg_ff(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

__device__ double g_ggg_ff(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

__device__ double g_gqq_ff(double zm, double zp, double rand) {
  return zm + (zp - zm) * rand;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SPLITTING FUNCTIONS BELOW - NO NEED TO EDIT
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

/**
 * splitting functions as a function - safer but less sophisticated
 * ----------------------------------------------------------------
 *
 * NB: HAPPY TO HEAR COMMENTS AND TIPS TO IMPROVE THIS COMPONENT
 *
 * this is a safer and more straightforward way to implement the splitting
 * functions for the shower. although the class-based approach is good for
 * CPU, in GPU many issues arise that mean that oop might not always be the
 * best strategy in coding. as a simpler approach, we will use switch-case
 * statements to select the correct splitting function.
 *
 * we have a lot of splitting functions:
 * - four CASES (FF, fi, if, ii)
 * - three or four possible TYPES (q->qg, q->gq, g->gg, g->qqbar, g->qbarq)
 * - five FLAVOURS of quarks (d, u, s, c, b) and each of their antiquarks
 * - at most, in total: 4 * (10 + 10 + 1 + 5) = 104 splitting functions
 *
 * so we need to organise ourselves with some kind of structure. as a first
 * attempt lets use four digit codes to identify the splitting functions:
 *
 * - 1st: CASE (FF, fi, if, ii) - 0, 1, 2, 3
 * - 2nd: SPLITTING (q->qg, q->gq, g->gg, g->qqbar, g->qbarq) - 0, 1, 2, 3, 4
 * - 3rd: is the emitter an ANTIPARTICLE - 0, 1 (quark, gluon is 0)
 * - 4th: FLAVOUR - 1, 2, 3, 4, 5; 0 for gluon
 *
 * examples:
 * - FF u -> u g = 0 0 0 2
 * - FF ubar -> ubar g = 0 0 1 2
 * - FF g -> g g = 0 2 0 0
 * - FF g -> c cbar = 0 3 0 4
 * - (Same as cbar c in FF, so divide sf value by 2)
 *
 * - fi u -> u g = 1 0 0 2
 * - fi ubar -> ubar g = 1 0 1 2
 * - fi g -> g g = 1 2 0 0
 * - fi g -> c cbar = 1 3 0 4
 * - (Same as cbar c in fi, so divide sf value by 2)
 *
 * - if d -> d g = 2 0 0 1
 * - if d -> g d = 2 1 0 1
 * - if sbar -> sbar g = 2 0 1 3
 * - if sbar -> g sbar = 2 1 1 3
 * - if g -> g g = 2 2 0 0
 * - if g -> u ubar = 2 3 0 2
 * - if g -> ubar u = 2 4 0 2
 *
 * - ii d -> d g = 3 0 0 1
 * - ii d -> g d = 3 1 0 1
 * - ii cbar -> cbar g = 3 0 1 4
 * - ii cbar -> g cbar = 3 1 1 4
 * - ii g -> g g = 3 2 0 0
 * - ii g -> b bbar = 3 3 0 5
 * - ii g -> bbar b = 3 4 0 5
 *
 * this way we can easily identify the splitting functions and select the
 * correct one using a switch-case statement. this can be used for value,
 * estimate, integral and generate_z functions.
 */

// -----------------------------------------------------------------------------
// Utility Functions to avoid magic numbers

// FF = 0, FI = 1, IF = 2, II = 3 - Not needed in GAPS 1
// __device__ int shower::get_splitting_case(int sf) const {
//   /**
//    * @brief Get the splitting case of the splitting function. This is the
//    first
//    * digit of the splitting function code
//    *
//    * @param sf: Splitting function code
//    * @return int: Splitting case
//    */

//   return sf / 1000;
// }

// q->qg = 0, q->gq = 1, g->gg = 2, g->qq = 3
__device__ int shower::get_splitting_type(int sf) const {
  /**
   * @brief Get the splitting type of the splitting function. This is the second
   * digit of the splitting function code
   *
   * @param sf: Splitting function code
   * @return int: Splitting type
   */

  return (sf % 1000) / 100;
}

// particle = 0, antiparticle = 1
__device__ int shower::is_emitter_antiparticle(int sf) const {
  /**
   * @brief Check if the emitter is an antiparticle. This is determined by the
   * third digit of the splitting function code
   *
   * @param sf: Splitting function code
   * @return int: 1 if antiparticle, 0 if particle
   */

  return (sf % 100) / 10;
}

// d = 1, u = 2, s = 3, c = 4, b = 5; gluon = 0
__device__ int shower::get_splitting_flavour(int sf) const {
  /**
   * @brief Get the flavour of the splitting function. This is the last digit of
   * the splitting function code
   *
   * @param sf: Splitting function code
   * @return int: Flavour of the splitting function
   */

  return sf % 10;
}

// -----------------------------------------------------------------------------
// Splitting Functions - Values

__device__ double shower::sf_value(double z, double y, int sf) const {
  /**
   * @brief Calculate the value of the splitting function for a given z and y
   *
   * @param z: z value of the splitting function
   * @param y: y value of the splitting function
   * @param sf: splitting function code
   * @return double: Splitting function value
   */

  switch (get_splitting_type(sf)) {
    case 0:
      return v_qqg_ff(z, y);
    case 2:
      return v_ggg_ff(z, y);
    case 3:
      return v_gqq_ff(z, y);
    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------
// Splitting Functions - Estimates

__device__ double shower::sf_estimate(double z, int sf) const {
  /**
   * @brief Calculate the integral of the splitting function for a given zm and
   * zp
   *
   * @param zm: lower z boundary of the integral
   * @param zp: upper z boundary of the integral
   * @param sf: splitting function code
   * @return double: Splitting function integral
   */

  switch (get_splitting_type(sf)) {
    case 0:
      return e_qqg_ff(z);
    case 2:
      return e_ggg_ff(z);
    case 3:
      return e_gqq_ff(z);
    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------
// Splitting Functions - Integrals

__device__ double shower::sf_integral(double zm, double zp, int sf) const {
  /**
   * @brief Generate a z value for the splitting function
   *
   * @param zm: lower z boundary of the integral
   * @param zp: upper z boundary of the integral
   * @param rand: random number between 0 and 1
   * @param sf: splitting function code
   * @return double: z value
   */

  switch (get_splitting_type(sf)) {
    case 0:
      return i_qqg_ff(zm, zp);
    case 2:
      return i_ggg_ff(zm, zp);
    case 3:
      return i_gqq_ff(zm, zp);
    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------
// Splitting Functions - Generate Z

__device__ double shower::sf_generate_z(double zm, double zp, double rand,
                                        int sf) const {
  /**
   * @brief Generate a z value for the splitting function
   *
   * @param zm: lower z boundary of the integral
   * @param zp: upper z boundary of the integral
   * @param rand: random number between 0 and 1
   * @param sf: splitting function code
   * @return double: z value
   */

  switch (get_splitting_type(sf)) {
    case 0:
      return g_qqg_ff(zm, zp, rand);
    case 2:
      return g_ggg_ff(zm, zp, rand);
    case 3:
      return g_gqq_ff(zm, zp, rand);
    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------

__device__ void shower::get_possible_splittings(int ij, int* splittings) const {
  /**
   * @brief Generate all possible splittings for a given dipole
   *
   * Why check so manY splitting functions, if 11 are enough?
   *
   * Need to determine the possible splittings rather than testing all
   * possible splittings. This is a major bottleneck in the code.
   *
   * Given a combination of emitter and spectator (FF, FI, IF, II), if the
   * emitter is a quark, there are only two possible splittings: q -> qg and
   * q -> gq. If the emitter is a gluon, there are six possible splittings:
   * g -> gg and ten g -> qq
   *
   * @param ij_pid: particle ID of the emitter
   * @param k_pid: particle ID of the spectator
   * @param ij_init: is the emitter an initial state particle
   * @param k_init: is the spectator an initial state particle
   * @param sf_codes: array to store the splitting function codes
   * @return int: number of possible splittings
   */

  // gluon splittings
  if (ij == 21) {
    splittings[0] = 200;
    splittings[1] = 301;
    splittings[2] = 302;
    splittings[3] = 303;
    splittings[4] = 304;
    splittings[5] = 305;
    // u, d, s only
    // splittings[4] = -1;
    // splittings[5] = -1;
  }
  // quark splittings
  else {
    if (ij < 0) {
      splittings[0] = -ij + 10;
      splittings[1] = -1;
    } else {
      splittings[0] = ij;
      splittings[1] = -1;
    }
  }
}

__device__ bool shower::validate_splitting(int ij, int sf, bool emt_init,
                                           bool spc_init) const {
  /**
   * @brief Validate the Splitting in the Select Winner Step
   *
   * Validate the Splitting in the Select Winner Step of Veto Algorithm
   * i.e. "is this splitting allowed for the given dipole?"
   *
   * UPDATE: Not used anymore, as we now use the generate_possible_splittings
   * to generate all possible splittings and then select the winner, and these
   * generated splittings are all valid by definition!
   *
   * @param ij: particle ID of the emitter
   * @param sf: splitting function code
   * @param emt_init: is the emitter an initial state particle
   * @param spc_init: is the spectator an initial state particle
   * @return bool: is the splitting valid
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II - Not needed in GAPS 1
  // int splitting_case = get_splitting_case(sf);

  // 0 = q->qg, 1 = q->gq, 2 = g->gg, 3 = g->qqbar, 4 = g->qbarq
  int splitting_type = get_splitting_type(sf);

  // 0 = particle, 1 = antiparticle
  int is_em_antip = is_emitter_antiparticle(sf);

  // 1 = d, 2 = u, 3 = s, 4 = c, 5 = b; 0 = gluon
  int splitting_flavour = get_splitting_flavour(sf);

  // skip if ij is a quark and the sf is not a quark sf (2nd digit), or
  // if ij is a gluon and the sf is not a gluon sf (2nd digit)
  if ((ij != 21 && splitting_type >= 2) || (ij == 21 && splitting_type < 2)) {
    return false;
  }

  // for quarks, skip if ij is a particle and sf is an antiparticle sf (3rd
  // digit), or if ij is an antiparticle and sf is a particle sf (3rd digit)
  if (ij != 21) {
    if ((ij < 0 && is_em_antip == 0) || (ij > 0 && is_em_antip == 1)) {
      return false;
    }
  }

  // for gluons, is_em_antip can never be 1
  else if (ij == 21 && is_em_antip == 1) {
    return false;
  }

  // q->qg case: Skip if the flavor of ij is different from the flavor of the
  // sf g->gg and g->qq case: No need to check the flavor
  if ((ij != 21 && abs(ij) != splitting_flavour)) {
    return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
// Convert Splitting Function to the Flavours Involved

__device__ void shower::sf_to_flavs(int sf, int* flavs) const {
  /**
   * @brief Convert the splitting function code to the flavours involved
   *
   * @param sf: splitting function code
   * @param flavs: array to store the flavours involved
   * @return int: number of flavours involved
   */

  // 0 = q->qg, 1 = q->gq, 2 = g->gg, 3 = g->qq, 4 = g->qbarq
  int splitting_type = get_splitting_type(sf);

  // 0 = particle, 1 = antiparticle
  int is_em_antip = is_emitter_antiparticle(sf);

  // 1 = d, 2 = u, 3 = s, 4 = c, 5 = b; 0 = gluon
  int splitting_flavour = get_splitting_flavour(sf);

  // q -> qg splittings ---------------------------------------

  if (splitting_type == 0) {
    if (is_em_antip == 0) {
      flavs[0] = splitting_flavour;
      flavs[1] = splitting_flavour;
      flavs[2] = 21;
    } else {
      flavs[0] = -1 * splitting_flavour;
      flavs[1] = -1 * splitting_flavour;
      flavs[2] = 21;
    }
  }

  // No q -> gq splittings in GAPS 1

  // g -> gg splittings ---------------------------------------

  if (splitting_type == 2) {
    flavs[0] = 21;
    flavs[1] = 21;
    flavs[2] = 21;
    return;
  }

  // g -> qq splittings ---------------------------------------

  if (splitting_type == 3) {
    flavs[0] = 21;
    flavs[1] = splitting_flavour;
    flavs[2] = -1 * splitting_flavour;
    return;
  }
}
