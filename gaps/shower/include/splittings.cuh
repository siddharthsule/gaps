#ifndef SPLITTINGS_CUH_
#define SPLITTINGS_CUH_

#include "qcd.cuh"

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

// Splitting Function Codes - Only FF for now
__constant__ int sfCodes[] = {0001, 0002, 0003, 0004, 0005, 0011, 0012, 0013,
                              0014, 0015, 0200, 0301, 0302, 0303, 0304, 0305};

__device__ double sfValue(double z, double y, int sfCode) {
  switch (sfCode) {
    // FF Splittings ---------------------------

    // FF q -> qg
    case 0001:
    case 0002:
    case 0003:
    case 0004:
    case 0005:
    case 0011:
    case 0012:
    case 0013:
    case 0014:
    case 0015:
      return kCF * (2.0 / (1.0 - z * (1.0 - y)) - (1.0 + z));
      break;

    // FF g -> gg
    case 0200:
      return kCA / 2.0 * (2.0 / (1.0 - z * (1.0 - y)) - 2.0 + z * (1.0 - z));
      break;

    // FF g -> qqbar
    case 0301:
    case 0302:
    case 0303:
    case 0304:
    case 0305:
      return kTR / 2.0 * (1.0 - 2.0 * z * (1.0 - z));
      break;
  }
}

__device__ double sfEstimate(double z, int sfCode) {
  switch (sfCode) {
    // FF Splittings ---------------------------

    // FF q -> qg
    case 0001:
    case 0002:
    case 0003:
    case 0004:
    case 0005:
    case 0011:
    case 0012:
    case 0013:
    case 0014:
    case 0015:
      return kCF * 2.0 / (1.0 - z);
      break;

    // FF g -> gg
    case 0200:
      return kCA / (1.0 - z);
      break;

    // FF g -> qqbar
    case 0301:
    case 0302:
    case 0303:
    case 0304:
    case 0305:
      return kTR / 2.0;
      break;
  }
}

__device__ double sfIntegral(double zm, double zp, int sfCode) {
  switch (sfCode) {
    // FF Splittings ---------------------------

    // FF q -> qg
    case 0001:
    case 0002:
    case 0003:
    case 0004:
    case 0005:
    case 0011:
    case 0012:
    case 0013:
    case 0014:
    case 0015:
      return kCF * 2.0 * log((1.0 - (1.0 - zp)) / (1.0 - zp));
      break;

    // FF g -> gg
    case 0200:
      return kCA * log((1.0 - (1.0 - zp)) / (1.0 - zp));
      break;

    // FF g -> qqbar
    case 0301:
    case 0302:
    case 0303:
    case 0304:
    case 0305:
      return kTR / 2.0 * (zp - (1.0 - zp));
      break;
  }
}

__device__ double sfGenerateZ(double zm, double zp, double rand, int sfCode) {
  switch (sfCode) {
    // FF Splittings ---------------------------

    // FF q -> qg
    case 0001:
    case 0002:
    case 0003:
    case 0004:
    case 0005:
    case 0011:
    case 0012:
    case 0013:
    case 0014:
    case 0015:
      return 1.0 + (zp - 1.0) * pow((1.0 - zm) / (1.0 - zp), rand);
      break;

    // FF g -> gg
    case 0200:
      return 1.0 + (zp - 1.0) * pow((1.0 - zm) / (1.0 - zp), rand);
      break;

    // FF g -> qqbar
    case 0301:
    case 0302:
    case 0303:
    case 0304:
    case 0305:
      return zm + (zp - zm) * rand;
      break;
  }
}

__device__ bool validateSplitting(int ij, int sf) {
  // Obtain the splitting function code
  int firstDigit = number / 1000;
  int secondDigit = (number / 100) % 10;
  int thirdDigit = (number / 10) % 10;
  int fourthDigit = number % 10;

  // Insert FF, FI, IF, II checks here
  // ---------------------------------

  // Skip if ij is a quark and the sf is not a quark sf (2nd digit), or
  // if ij is a gluon and the sf is not a gluon sf (2nd digit)
  if ((ij != 21 && secondDigit >= 2) || (ij == 21 && secondDigit < 2)) {
    return false;
  }

  // Skip if ij is a particle and sf is an antiparticle sf (3rd digit), or
  // if ij is an antiparticle and sf is a particle sf (3rd digit)
  if ((ij < 0 && thirdDigit == 0) || (ij > 0 && sf thirdDigit == 1)) {
    return false;
  }

  // Skip if the flavor of ij is different from the flavor of the sf
  if ((ij != 21 && abs(ij) != fourthDigit) || (ij == 21 && fourthDigit != 0)) {
    return false;
  }

  return true;
}

#endif  // SPLITTINGS_CUH_
