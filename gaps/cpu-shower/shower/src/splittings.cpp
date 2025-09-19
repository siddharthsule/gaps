#include "shower.h"

/**
 * Here you can define and edit the splitting functions for the shower. Make
 * sure to update the estimates, integrals and generate_z functions as well,
 * or else the vetoalgorithm will not work correctly.
 *
 * the functions are defined as (type)_(splitting)_(case):
 * - Type: V = Value, E = Estimate, I = Integral, G = Generate Z
 * - Splitting: q -> qg, q -> gq, g -> gg, g -> qq
 * - Case: FF, FI, IF, II
 *
 * To understand how these are Implemented, see the code after the functions.
 */

// -----------------------------------------------------------------------------
// FF splittings: z and y

// Values --------------------------------

double v_qqg_ff(double z, double y) {
  return k_cf * (2. / (1. - z * (1. - y)) - (1. + z));
}

double v_ggg_ff(double z, double y) {
  return k_ca / 2. * (2. / (1. - z * (1. - y)) - 2. + z * (1. - z));
}

double v_gqq_ff(double z, double y) {
  return k_tr / 2. * (1. - 2. * z * (1. - z));
}

// Estimates ------------------------------

double e_qqg_ff(double z) { return k_cf * 2. / (1. - z); }

double e_ggg_ff(double z) { return k_ca / (1. - z); }

double e_gqq_ff(double z) { return k_tr / 2.; }

// Integrals ------------------------------

double i_qqg_ff(double zm, double zp) {
  return k_cf * 2. * log((1. - zm) / (1. - zp));
}

double i_ggg_ff(double zm, double zp) {
  return k_ca * log((1. - zm) / (1. - zp));
}

double i_gqq_ff(double zm, double zp) { return k_tr / 2. * (zp - zm); }

// Generate Z ------------------------------

double g_qqg_ff(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

double g_ggg_ff(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

double g_gqq_ff(double zm, double zp, double rand) {
  return zm + (zp - zm) * rand;
}

// -----------------------------------------------------------------------------
// FI splittings: z and x (here, y)

// Values --------------------------------

double v_qqg_fi(double z, double y) {
  return k_cf * (2. / (1. - z + (1. - y)) - (1. + z));
}

double v_ggg_fi(double z, double y) {
  return k_ca / 2. * (2. / (1. - z + (1. - y)) - 2. + z * (1. - z));
}

double v_gqq_fi(double z, double y) {
  return k_tr / 2. * (1. - 2. * z * (1. - z));
}

// Estimates ------------------------------

double e_qqg_fi(double z) { return k_cf * 2. / (1. - z); }

double e_ggg_fi(double z) { return k_ca / (1. - z); }

double e_gqq_fi(double z) { return k_tr / 2.; }

// Integrals ------------------------------

double i_qqg_fi(double zm, double zp) {
  return k_cf * 2. * log((1. - zm) / (1. - zp));
}

double i_ggg_fi(double zm, double zp) {
  return k_ca * log((1. - zm) / (1. - zp));
}

double i_gqq_fi(double zm, double zp) { return k_tr / 2. * (zp - zm); }

// Generate Z ------------------------------

double g_qqg_fi(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

double g_ggg_fi(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

double g_gqq_fi(double zm, double zp, double rand) {
  return zm + (zp - zm) * rand;
}

// -----------------------------------------------------------------------------
// IF splittings: x (here, z) and u (here, y)

// Values --------------------------------

double v_qqg_if(double z, double y) {
  return k_cf * (2. / (1. - z + y) - (1. + z));
}

// q -> gq backwards is g -> q qbar!
double v_qgq_if(double z, double y) { return k_tr * (1. - 2. * z * (1. - z)); }

double v_ggg_if(double z, double y) {
  // return k_ca * (1. / (1. - z + y) + (1. - z) / z - 1. + z * (1. - z));
  return k_ca * (1. / (1. - z + y) + 1. / z - 2. + z * (1. - z));
}

// g -> q q backwards is q -> gq!
// Also now g -> q qbar != g -> qbar q
double v_gqqbar_if(double z, double y) {
  return k_cf * (z + 2. * (1. - z) / z);
}

double v_gqbarq_if(double z, double y) {
  return k_cf * (z + 2. * (1. - z) / z);
}

// Estimates ------------------------------

double e_qqg_if(double z) { return k_cf * 2. / (1. - z); }

double e_qgq_if(double z) { return k_tr; }

double e_ggg_if(double z) { return k_ca * (1 / (1. - z) + 1 / z); }

double e_gqqbar_if(double z) { return k_cf * 2. / z; }

double e_gqbarq_if(double z) { return k_cf * 2. / z; }

// Integrals ------------------------------

double i_qqg_if(double zm, double zp) {
  return k_cf * 2. * log((1. - zm) / (1. - zp));
}

double i_qgq_if(double zm, double zp) { return k_tr * (zp - zm); }

double i_ggg_if(double zm, double zp) {
  return k_ca * log((zp * (1. - zm)) / (zm * (1. - zp)));
}

double i_gqqbar_if(double zm, double zp) { return k_cf * 2 * log(zp / zm); }

double i_gqbarq_if(double zm, double zp) { return k_cf * 2 * log(zp / zm); }

// Generate Z ------------------------------

double g_qqg_if(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

double g_qgq_if(double zm, double zp, double rand) {
  return zm + (zp - zm) * rand;
}

double g_ggg_if(double zm, double zp, double rand) {
  double m = pow(zm / (1.0 - zm), rand - 1) * pow(zp / (1.0 - zp), -rand);
  return 1. / (1. + m);
}

double g_gqqbar_if(double zm, double zp, double rand) {
  return zm * pow(zp / zm, rand);
}

double g_gqbarq_if(double zm, double zp, double rand) {
  return zm * pow(zp / zm, rand);
}

// -----------------------------------------------------------------------------
// II splittings: x (here, z) and v (here, y)

// Values --------------------------------

double v_qqg_ii(double z, double y) {
  return k_cf * (2. / (1. - z) - (1. + z));
}

// q -> gq backwards is g -> q qbar!
double v_qgq_ii(double z, double y) { return k_tr * (1. - 2. * z * (1. - z)); }

double v_ggg_ii(double z, double y) {
  // return k_ca * (1. / (1. - z) + (1. - z) / z - 1. + z * (1. - z));
  return k_ca * (1. / (1. - z) + 1. / z - 2. + z * (1. - z));
}

// g -> q q backwards is q -> gq!
// Also now g -> q qbar != g -> qbar q
double v_gqqbar_ii(double z, double y) {
  return k_cf / 2. * (z + 2. * (1. - z) / z);
}

double v_gqbarq_ii(double z, double y) {
  return k_cf / 2. * (z + 2. * (1. - z) / z);
}

// Estimates ------------------------------

double e_qqg_ii(double z) { return k_cf * 2. / (1. - z); }

double e_qgq_ii(double z) { return k_tr; }

double e_ggg_ii(double z) { return k_ca * (1 / (1. - z) + 1 / z); }

double e_gqqbar_ii(double z) { return k_cf * 2. / z; }

double e_gqbarq_ii(double z) { return k_cf * 2. / z; }

// Integrals ------------------------------

double i_qqg_ii(double zm, double zp) {
  return k_cf * 2. * log((1. - zm) / (1. - zp));
}

double i_qgq_ii(double zm, double zp) { return k_tr * (zp - zm); }

double i_ggg_ii(double zm, double zp) {
  return k_ca * log((zp * (1. - zm)) / (zm * (1. - zp)));
}

double i_gqqbar_ii(double zm, double zp) { return k_cf * 2. * log(zp / zm); }

double i_gqbarq_ii(double zm, double zp) { return k_cf * 2. * log(zp / zm); }

// Generate Z ------------------------------

double g_qqg_ii(double zm, double zp, double rand) {
  return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
}

double g_qgq_ii(double zm, double zp, double rand) {
  return zm + (zp - zm) * rand;
}

double g_ggg_ii(double zm, double zp, double rand) {
  double m = pow(zm / (1.0 - zm), rand - 1) * pow(zp / (1.0 - zp), -rand);
  return 1. / (1. + m);
}

double g_gqqbar_ii(double zm, double zp, double rand) {
  return zm * pow(zp / zm, rand);
}

double g_gqbarq_ii(double zm, double zp, double rand) {
  return zm * pow(zp / zm, rand);
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

// FF = 0, FI = 1, IF = 2, II = 3
int shower::get_splitting_case(int sf) const {
  /**
   * @brief Get the splitting case of the splitting function. This is the first
   * digit of the splitting function code
   *
   * @param sf: Splitting function code
   * @return int: Splitting case
   */

  return sf / 1000;
}

bool shower::is_ff(int sf) const {
  /**
   * @brief Check if the splitting function is a FF splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if FF, False otherwise
   */

  return sf / 1000 == 0;
}

bool shower::is_fi(int sf) const {
  /**
   * @brief Check if the splitting function is a FI splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if FI, False otherwise
   */

  return sf / 1000 == 1;
}

bool shower::is_if(int sf) const {
  /**
   * @brief Check if the splitting function is an IF splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if IF, False otherwise
   */

  return sf / 1000 == 2;
}

bool shower::is_ii(int sf) const {
  /**
   * @brief Check if the splitting function is an II splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if II, False otherwise
   */

  return sf / 1000 == 3;
}
bool shower::is_fsr(int sf) const {
  /**
   * @brief Check if the splitting function is an FSR splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if FSR, False otherwise
   */

  return sf / 1000 <= 1;
}

bool shower::is_isr(int sf) const {
  /**
   * @brief Check if the splitting function is an ISR splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if ISR, False otherwise
   */

  return sf / 1000 > 1;
}

// q->qg = 0, q->gq = 1, g->gg = 2, g->qq = 3
int shower::get_splitting_type(int sf) const {
  /**
   * @brief Get the splitting type of the splitting function. This is the second
   * digit of the splitting function code
   *
   * @param sf: Splitting function code
   * @return int: Splitting type
   */

  return (sf % 1000) / 100;
}

bool shower::is_q2qg(int sf) const {
  /**
   * @brief Check if the splitting function is a q->qg splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if q->qg, False otherwise
   */

  return (sf % 1000) / 100 == 0;
}

bool shower::is_q2gq(int sf) const {
  /**
   * @brief Check if the splitting function is a q->gq splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if q->gq, False otherwise
   */

  return (sf % 1000) / 100 == 1;
}

bool shower::is_g2gg(int sf) const {
  /**
   * @brief Check if the splitting function is a g->gg splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if g->gg, False otherwise
   */

  return (sf % 1000) / 100 == 2;
}

bool shower::is_g2qqbar(int sf) const {
  /**
   * @brief Check if the splitting function is a g->qq splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if g->qq, False otherwise
   */

  return (sf % 1000) / 100 == 3;
}

bool shower::is_g2qbarq(int sf) const {
  /**
   * @brief Check if the splitting function is a q->gqbar splitting
   *
   * @param sf: Splitting function code
   * @return bool: True if q->gqbar, False otherwise
   */

  return (sf % 1000) / 100 == 4;
}

// particle = 0, antiparticle = 1
int shower::is_emitter_antiparticle(int sf) const {
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
int shower::get_splitting_flavour(int sf) const {
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

double shower::sf_value(double z, double y, int sf) const {
  /**
   * @brief Calculate the value of the splitting function for a given z and y
   *
   * @param z: z value of the splitting function
   * @param y: y value of the splitting function
   * @param sf: splitting function code
   * @return double: Splitting function value
   */

  switch (get_splitting_case(sf)) {
    // FF
    case 0:
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

    // FI
    case 1:
      switch (get_splitting_type(sf)) {
        case 0:
          return v_qqg_fi(z, y);
        case 2:
          return v_ggg_fi(z, y);
        case 3:
          return v_gqq_fi(z, y);
        default:
          return 0;
      }

    // IF
    case 2:
      switch (get_splitting_type(sf)) {
        case 0:
          return v_qqg_if(z, y);
        case 1:
          return v_qgq_if(z, y);
        case 2:
          return v_ggg_if(z, y);
        case 3:
          return v_gqqbar_if(z, y);
        case 4:
          return v_gqbarq_if(z, y);
        default:
          return 0;
      }

    // II
    case 3:
      switch (get_splitting_type(sf)) {
        case 0:
          return v_qqg_ii(z, y);
        case 1:
          return v_qgq_ii(z, y);
        case 2:
          return v_ggg_ii(z, y);
        case 3:
          return v_gqqbar_ii(z, y);
        case 4:
          return v_gqbarq_ii(z, y);
        default:
          return 0;
      }

    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------
// Splitting Functions - Estimates

double shower::sf_estimate(double z, int sf) const {
  /**
   * @brief Calculate the estimate of the splitting function for a given z
   *
   * @param z: z value of the splitting function
   * @param sf: splitting function code
   * @return double: Splitting function estimate
   */

  switch (get_splitting_case(sf)) {
    // FF
    case 0:
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

    // FI
    case 1:
      switch (get_splitting_type(sf)) {
        case 0:
          return e_qqg_fi(z);
        case 2:
          return e_ggg_fi(z);
        case 3:
          return e_gqq_fi(z);
        default:
          return 0;
      }

    // IF
    case 2:
      switch (get_splitting_type(sf)) {
        case 0:
          return e_qqg_if(z);
        case 1:
          return e_qgq_if(z);
        case 2:
          return e_ggg_if(z);
        case 3:
          return e_gqqbar_if(z);
        case 4:
          return e_gqbarq_if(z);
        default:
          return 0;
      }

    // II
    case 3:
      switch (get_splitting_type(sf)) {
        case 0:
          return e_qqg_ii(z);
        case 1:
          return e_qgq_ii(z);
        case 2:
          return e_ggg_ii(z);
        case 3:
          return e_gqqbar_ii(z);
        case 4:
          return e_gqbarq_ii(z);
        default:
          return 0;
      }

    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------
// Splitting Functions - Integrals

double shower::sf_integral(double zm, double zp, int sf) const {
  /**
   * @brief Calculate the integral of the splitting function for a given zm and
   * zp
   *
   * @param zm: lower z boundary of the integral
   * @param zp: upper z boundary of the integral
   * @param sf: splitting function code
   * @return double: Splitting function integral
   */

  switch (get_splitting_case(sf)) {
    // FF
    case 0:
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

    // FI
    case 1:
      switch (get_splitting_type(sf)) {
        case 0:
          return i_qqg_fi(zm, zp);
        case 2:
          return i_ggg_fi(zm, zp);
        case 3:
          return i_gqq_fi(zm, zp);
        default:
          return 0;
      }

    // IF
    case 2:
      switch (get_splitting_type(sf)) {
        case 0:
          return i_qqg_if(zm, zp);
        case 1:
          return i_qgq_if(zm, zp);
        case 2:
          return i_ggg_if(zm, zp);
        case 3:
          return i_gqqbar_if(zm, zp);
        case 4:
          return i_gqbarq_if(zm, zp);
        default:
          return 0;
      }

    // II
    case 3:
      switch (get_splitting_type(sf)) {
        case 0:
          return i_qqg_ii(zm, zp);
        case 1:
          return i_qgq_ii(zm, zp);
        case 2:
          return i_ggg_ii(zm, zp);
        case 3:
          return i_gqqbar_ii(zm, zp);
        case 4:
          return i_gqbarq_ii(zm, zp);

        default:
          return 0;
      }

    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------
// Splitting Functions - Generate Z

double shower::sf_generate_z(double zm, double zp, double rand, int sf) const {
  /**
   * @brief Generate a z value for the splitting function
   *
   * @param zm: lower z boundary of the integral
   * @param zp: upper z boundary of the integral
   * @param rand: random number between 0 and 1
   * @param sf: splitting function code
   * @return double: z value
   */

  switch (get_splitting_case(sf)) {
    // FF
    case 0:
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

    // FI
    case 1:
      switch (get_splitting_type(sf)) {
        case 0:
          return g_qqg_fi(zm, zp, rand);
        case 2:
          return g_ggg_fi(zm, zp, rand);
        case 3:
          return g_gqq_fi(zm, zp, rand);
        default:
          return 0;
      }

    // IF
    case 2:
      switch (get_splitting_type(sf)) {
        case 0:
          return g_qqg_if(zm, zp, rand);
        case 1:
          return g_qgq_if(zm, zp, rand);
        case 2:
          return g_ggg_if(zm, zp, rand);
        case 3:
          return g_gqqbar_if(zm, zp, rand);
        case 4:
          return g_gqbarq_if(zm, zp, rand);
        default:
          return 0;
      }

    // II
    case 3:
      switch (get_splitting_type(sf)) {
        case 0:
          return g_qqg_ii(zm, zp, rand);
        case 1:
          return g_qgq_ii(zm, zp, rand);
        case 2:
          return g_ggg_ii(zm, zp, rand);
        case 3:
          return g_gqqbar_ii(zm, zp, rand);
        case 4:
          return g_gqbarq_ii(zm, zp, rand);
        default:
          return 0;
      }

    default:
      return 0;
  }
}

// -----------------------------------------------------------------------------

void shower::generate_possible_splittings(int ij_pid, int k_pid, bool ij_init,
                                          bool k_init, int* sf_codes) const {
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

  // Initialise the array with null splittings
  for (int i = 0; i < 11; i++) {
    sf_codes[i] = -1;
  }

  // No need to fill all elements
  int possible_splittings;
  if (ij_pid == 21) {
    if (ij_init) {
      possible_splittings = 11;
    } else {
      possible_splittings = 6;
    }
  } else {
    if (ij_init) {
      possible_splittings = 2;
    } else {
      possible_splittings = 1;
    }
  }

  // Sort out the splitting case (1st) code
  for (int i = 0; i < possible_splittings; i++) {
    // FF
    if (!ij_init && !k_init) {
      sf_codes[i] = 0;
    }
    // FI
    if (!ij_init && k_init) {
      sf_codes[i] = 1000;
    }
    // IF
    if (ij_init && !k_init) {
      sf_codes[i] = 2000;
    }
    // II
    if (ij_init && k_init) {
      sf_codes[i] = 3000;
    }
  }

  // Sort out the type (2nd), antiparticle (3rd) and flavour (4th) codes
  if (ij_pid == 21) {
    sf_codes[0] += 200;  // g -> gg
    sf_codes[1] += 301;  // g -> d dbar
    sf_codes[2] += 302;  // g -> u ubar
    sf_codes[3] += 303;  // g -> s sbar
    sf_codes[4] += 304;  // g -> c cbar
    sf_codes[5] += 305;  // g -> b bbar
    // Add g -> qbar q if initial gluon
    if (ij_init) {
      sf_codes[6] += 401;   // g -> ubar u
      sf_codes[7] += 402;   // g -> dbar d
      sf_codes[8] += 403;   // g -> sbar s
      sf_codes[9] += 404;   // g -> cbar c
      sf_codes[10] += 405;  // g -> bbar b
    }
  }

  // Add q -> gq (and q -> qg if initial quark)
  else {
    if (ij_pid > 0) {
      sf_codes[0] += ij_pid;  // q -> qg
      if (ij_init) {
        sf_codes[1] += ij_pid + 100;  // q -> gq
      }
    } else {
      sf_codes[0] += -ij_pid + 10;  // q -> qg
      if (ij_init) {
        sf_codes[1] += -ij_pid + 10 + 100;  // q -> gq
      }
    }
  }
}

bool shower::validate_splitting(int ij, int sf, bool emt_init,
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

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  // 0 = q->qg, 1 = q->gq, 2 = g->gg, 3 = g->qqbar, 4 = g->qbarq
  int splitting_type = get_splitting_type(sf);

  // 0 = particle, 1 = antiparticle
  int is_em_antip = is_emitter_antiparticle(sf);

  // 1 = d, 2 = u, 3 = s, 4 = c, 5 = b; 0 = gluon
  int splitting_flavour = get_splitting_flavour(sf);

  // The ! means that it's only correct if all three conditions are met
  // FF
  if (splitting_case == 0) {
    if (!(!emt_init && !spc_init)) {
      return false;
    }
  }
  // FI
  else if (splitting_case == 1) {
    if (!(!emt_init && spc_init)) {
      return false;
    }
  }
  // IF
  else if (splitting_case == 2) {
    if (!(emt_init && !spc_init)) {
      return false;
    }
  }
  // II
  else if (splitting_case == 3) {
    if (!(emt_init && spc_init)) {
      return false;
    }
  }

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

void shower::sf_to_flavs(int sf, int* flavs) const {
  /**
   * @brief Convert the splitting function code to the flavours involved
   *
   * @param sf: splitting function code
   * @param flavs: array to store the flavours involved
   * @return int: number of flavours involved
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

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

  // q -> gq splittings ---------------------------------------

  // Only occurs if first_digit > 1 (IF, II)
  // q -> gq backwards is g -> q qbar!
  if (splitting_type == 1) {
    if (is_em_antip == 0) {
      flavs[0] = splitting_flavour;
      flavs[1] = 21;
      flavs[2] = -1 * splitting_flavour;
    } else {
      flavs[0] = -1 * splitting_flavour;
      flavs[1] = 21;
      flavs[2] = splitting_flavour;
    }
  }

  // g -> gg splittings ---------------------------------------

  if (splitting_type == 2) {
    flavs[0] = 21;
    flavs[1] = 21;
    flavs[2] = 21;
    return;
  }

  // g -> qq splittings ---------------------------------------

  if (splitting_type == 3) {
    // FF, FI
    if (splitting_case < 2) {
      flavs[0] = 21;
      flavs[1] = splitting_flavour;
      flavs[2] = -1 * splitting_flavour;
      return;
    }
    // IF, II
    // g -> q qbar backwards is qbar -> g qbar!
    else {
      flavs[0] = 21;
      flavs[1] = -1 * splitting_flavour;
      flavs[2] = -1 * splitting_flavour;
      return;
    }
  }

  if (splitting_type == 4) {
    // IF, II
    // g -> qbar q backwards is q -> g q!
    flavs[0] = 21;
    flavs[1] = splitting_flavour;
    flavs[2] = splitting_flavour;
  }
}

// -----------------------------------------------------------------------------
// Convert sf code to Text for easier debugging

void shower::sf_to_text(int sf, char* text) const {
  /**
   * @brief Convert the splitting function code to text for easy debugging
   *
   * @param sf: splitting function code
   * @param text: array to store the text
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = sf / 1000;

  // 0 = q->qg, 1 = q->gq, 2 = g->gg, 3 = g->qq
  int splitting_type = (sf % 1000) / 100;

  // Convert Case
  switch (splitting_case) {
    case 0:
      text[0] = 'F';
      text[1] = 'F';
      break;
    case 1:
      text[0] = 'F';
      text[1] = 'I';
      break;
    case 2:
      text[0] = 'I';
      text[1] = 'F';
      break;
    case 3:
      text[0] = 'I';
      text[1] = 'I';
      break;
    default:
      text[0] = 'E';
      text[1] = 'S';
      break;
  }

  // Extra space
  text[2] = ' ';

  // Convert Type
  switch (splitting_type) {
    case 0:
      text[3] = 'q';
      text[4] = '-';
      text[5] = '>';
      text[6] = 'q';
      text[7] = 'g';
      break;
    case 1:
      text[3] = 'q';
      text[4] = '-';
      text[5] = '>';
      text[6] = 'g';
      text[7] = 'q';
      break;
    case 2:
      text[3] = 'g';
      text[4] = '-';
      text[5] = '>';
      text[6] = 'g';
      text[7] = 'g';
      break;
    case 3:
      text[3] = 'g';
      text[4] = '-';
      text[5] = '>';
      text[6] = 'q';
      text[7] = 'q';
      break;
    default:
      text[3] = 'E';
      text[4] = 'S';
      text[5] = ' ';
      text[6] = ' ';
      text[7] = ' ';
      break;
  }

  // Extra space
  text[8] = ' ';

  // Extra character - ignored
  text[9] = '\0';
}