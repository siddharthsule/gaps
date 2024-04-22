#include "shower.h"

/**
 * Splitting Functions as a function - safer but less sophisticated
 * ----------------------------------------------------------------
 *
 * This is a safer and more straightforward way to implement the splitting
 * functions for the shower. Although the class-based approach is good for
 * C++, in CUDA many issues arise that mean that OOP might not always be the
 * best strategy in coding. As a simpler approach, we will use switch-case
 * statements to select the correct splitting function.
 *
 * We have a LOT of splitting functions:
 * - Four types (FF, FI, IF, II)
 * - Three or Four Possible DGLAP Splittings (q->qg, q->gq, g->gg, g->qq)
 * - Five Flavours of Quarks (d, u, s, c, b) and each of their antiquarks
 * - At most, In total: 4 * (10 + 10 + 1 + 5) = 104 splitting functions
 *
 * So we need to organise ourselves with some kind of structure. As a first
 * attempt lets use four digit codes to identify the splitting functions:
 *
 * - 1st digit: Type of Split-Spect (FF, FI, IF, II) - 0, 1, 2, 3
 * - 2nd digit: Type of DGLAP (q->qg, q->gq, g->gg, g->qq) - 0, 1, 2, 3
 * - 3rd digit: Emitter is a Particle or Antiparticle - 0, 1 (gluon is 0)
 * - 4th digit: Flavor of the Emitter - 1, 2, 3, 4, 5; 0 for gluon
 *
 * Examples:
 * - FF u -> ug = 0 0 0 2
 * - FF ubar -> ubar g = 0 0 1 2
 * - FF g -> gg = 0 2 0 0
 * - FF g -> ccbar = 0 3 0 4
 *
 * - FI u -> ug = 1 0 0 2
 * - FI g -> ccbar = 1 3 0 4
 *
 * - IF d -> dg = 2 0 0 1
 * - IF d -> gd = 2 1 0 1
 * - IF sbar -> sbar g = 2 0 1 3
 * - IF g -> uubar = 2 3 0 2
 *
 * - II g -> gg = 3 2 0 0
 * - II g -> bbbar = 3 3 0 5
 *
 * This way we can easily identify the splitting functions and select the
 * correct one using a switch-case statement. This can be used for value,
 * estimate, integral and generateZ functions.
 */

double Shower::sfValue(double z, double y, int sf) {
  switch (sf) {
    // FF Splittings ---------------------------

    // FF q -> qg
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
      return kCF * (2. / (1. - z * (1. - y)) - (1. + z));
      break;

    // FF g -> gg
    case 200:
      return kCA / 2. * (2. / (1. - z * (1. - y)) - 2. + z * (1. - z));
      break;

    // FF g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return kTR / 2. * (1. - 2. * z * (1. - z));
      break;
  }
  return 0.;
}

double Shower::sfEstimate(double z, int sf) {
  switch (sf) {
    // FF Splittings ---------------------------

    // FF q -> qg
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
      return kCF * 2. / (1. - z);
      break;

    // FF g -> gg
    case 200:
      return kCA / (1. - z);
      break;

    // FF g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return kTR / 2.;
      break;
  }
  return 0.;
}

double Shower::sfIntegral(double zm, double zp, int sf) {
  switch (sf) {
    // FF Splittings ---------------------------

    // FF q -> qg
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
      return kCF * 2. * log((1. - zm) / (1. - zp));
      break;

    // FF g -> gg
    case 200:
      return kCA * log((1. - zm) / (1. - zp));
      break;

    // FF g -> qqbar
    case 301:
    case 302:
    case 303:
    case 304:
    case 305:
      return kTR / 2. * (zp - zm);
      break;
  }
  return 0.;
}

double Shower::sfGenerateZ(double zm, double zp, double rand, int sf) {
  switch (sf) {
    // FF Splittings ---------------------------

    // FF q -> qg
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

    // FF g -> gg
    case 200:
      return 1. + (zp - 1.) * pow((1. - zm) / (1. - zp), rand);
      break;

    // FF g -> qqbar
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

bool Shower::validateSplitting(int ij, int sf) {
  // Obtain the splitting function code
  // int firstDigit = sf / 1000;
  int secondDigit = (sf / 100) % 10;
  int thirdDigit = (sf / 10) % 10;
  int fourthDigit = sf % 10;

  // Insert FF, FI, IF, II checks here
  // ---------------------------------

  // Skip if ij is a quark and the sf is not a quark sf (2nd digit), or
  // if ij is a gluon and the sf is not a gluon sf (2nd digit)
  if ((ij != 21 && secondDigit >= 2) || (ij == 21 && secondDigit < 2)) {
    return false;
  }

  // Skip if ij is a particle and sf is an antiparticle sf (3rd digit), or
  // if ij is an antiparticle and sf is a particle sf (3rd digit)
  if ((ij < 0 && thirdDigit == 0) || (ij > 0 && thirdDigit == 1)) {
    return false;
  }

  // Skip if the flavor of ij is different from the flavor of the sf
  if ((ij != 21 && abs(ij) != fourthDigit) || (ij == 21 && fourthDigit != 0)) {
    return false;
  }

  return true;
}

void Shower::sfToFlavs(int sf, int* flavs) {
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
