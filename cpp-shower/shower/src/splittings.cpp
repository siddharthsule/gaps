#include "shower.h"

/**
 * splitting functions as a function - safer but less sophisticated
 * ----------------------------------------------------------------
 *
 * this is a safer and more straightforward way to implement the splitting
 * functions for the shower. although the class-based approach is good for
 * c++, in cuda many issues arise that mean that oop might not always be the
 * best strategy in coding. as a simpler approach, we will use switch-case
 * statements to select the correct splitting function.
 *
 * we have a lot of splitting functions:
 * - four types (ff, fi, if, ii)
 * - three or four possible dglap splittings (q->qg, q->gq, g->gg, g->qq)
 * - five flavours of quarks (d, u, s, c, b) and each of their antiquarks
 * - at most, in total: 4 * (10 + 10 + 1 + 5) = 104 splitting functions
 *
 * so we need to organise ourselves with some kind of structure. as a first
 * attempt lets use four digit codes to identify the splitting functions:
 *
 * - 1st digit: type of split-spect (ff, fi, if, ii) - 0, 1, 2, 3
 * - 2nd digit: type of dglap (q->qg, q->gq, g->gg, g->qq) - 0, 1, 2, 3
 * - 3rd digit: emitter is a particle or antiparticle - 0, 1 (gluon is 0)
 * - 4th digit: flavor of the emitter - 1, 2, 3, 4, 5; 0 for gluon
 *
 * examples:
 * - ff u -> ug = 0 0 0 2
 * - ff ubar -> ubar g = 0 0 1 2
 * - ff g -> gg = 0 2 0 0
 * - ff g -> ccbar = 0 3 0 4
 *
 * - fi u -> ug = 1 0 0 2
 * - fi g -> ccbar = 1 3 0 4
 *
 * - if d -> dg = 2 0 0 1
 * - if d -> gd = 2 1 0 1
 * - if sbar -> sbar g = 2 0 1 3
 * - if g -> uubar = 2 3 0 2
 *
 * - ii g -> gg = 3 2 0 0
 * - ii g -> bbbar = 3 3 0 5
 *
 * this way we can easily identify the splitting functions and select the
 * correct one using a switch-case statement. this can be used for value,
 * estimate, integral and generate_z functions.
 */

double shower::sf_value(double z, double y, int sf) {
  switch (sf) {
    // ff splittings ---------------------------

    // ff q -> qg
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
      return k_cf * (2. / (1. - z * (1. - y)) - (1. + z));
      break;

    // ff g -> gg
    case 200:
      return k_ca / 2. * (2. / (1. - z * (1. - y)) - 2. + z * (1. - z));
      break;

    // ff g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return k_tr / 2. * (1. - 2. * z * (1. - z));
      break;
  }
  return 0.;
}

double shower::sf_estimate(double z, int sf) {
  switch (sf) {
    // ff splittings ---------------------------

    // ff q -> qg
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
      return k_cf * 2. / (1. - z);
      break;

    // ff g -> gg
    case 200:
      return k_ca / (1. - z);
      break;

    // ff g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return k_tr / 2.;
      break;
  }
  return 0.;
}

double shower::sf_integral(double zm, double zp, int sf) {
  switch (sf) {
    // ff splittings ---------------------------

    // ff q -> qg
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
      return k_cf * 2. * log((1. - zm) / (1. - zp));
      break;

    // ff g -> gg
    case 200:
      return k_ca * log((1. - zm) / (1. - zp));
      break;

    // ff g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return k_tr / 2. * (zp - zm);
      break;
  }
  return 0.;
}

double shower::sf_generate_z(double zm, double zp, double rand, int sf) {
  switch (sf) {
    // ff splittings ---------------------------

    // ff q -> qg
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
      return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
      break;

    // ff g -> gg
    case 200:
      return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
      break;

    // ff g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return zm + (zp - zm) * rand;
      break;
  }
  return 0.;
}

bool shower::validate_splitting(int ij, int sf) {
  // obtain the splitting function code
  // int first_digit = sf / 1000;
  int second_digit = (sf / 100) % 10;
  int third_digit = (sf / 10) % 10;
  int fourth_digit = sf % 10;

  // insert ff, fi, if, ii checks here
  // ---------------------------------

  // skip if ij is a quark and the sf is not a quark sf (2nd digit), or
  // if ij is a gluon and the sf is not a gluon sf (2nd digit)
  if ((ij != 21 && second_digit >= 2) || (ij == 21 && second_digit < 2)) {
    return false;
  }

  // skip if ij is a particle and sf is an antiparticle sf (3rd digit), or
  // if ij is an antiparticle and sf is a particle sf (3rd digit)
  if ((ij < 0 && third_digit == 0) || (ij > 0 && third_digit == 1)) {
    return false;
  }

  // skip if the flavor of ij is different from the flavor of the sf
  if ((ij != 21 && abs(ij) != fourth_digit) ||
      (ij == 21 && fourth_digit != 0)) {
    return false;
  }

  return true;
}

void shower::sf_to_flavs(int sf, int* flavs) {
  if (sf < 16) {
    if (sf < 6) {
      flavs[0] = sf;
      flavs[1] = sf;
      flavs[2] = 21;
    } else {
      flavs[0] = -1 * (sf - 10);
      flavs[1] = -1 * (sf - 10);
      flavs[2] = 21;
    }
  } else if (sf == 200) {
    flavs[0] = 21;
    flavs[1] = 21;
    flavs[2] = 21;
  } else if (sf < 306) {
    flavs[0] = 21;
    flavs[1] = sf - 300;
    flavs[2] = -1 * (sf - 300);
  }
}
