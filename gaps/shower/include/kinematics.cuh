#ifndef SHOWER_KINEMATICS_CUH_
#define SHOWER_KINEMATICS_CUH_

#include "qcd.cuh"

/**
 * Why is this function in a header file? See Vec4.cuh
 */

__device__ void MakeKinematics(Vec4 *kinematics, const double z, const double y,
                               const double phi, const Vec4 pijt,
                               const Vec4 pkt) {
  Vec4 Q = pijt + pkt;

  // Generating the Momentum (0, kt1, kt2, 0)
  double rkt = sqrt(Q.M2() * y * z * (1. - z));

  Vec4 kt1 = pijt.Cross(pkt);
  if (kt1.P() < 1.e-6) {
    Vec4 xaxis(0., 1., 0., 0.);
    kt1 = pijt.Cross(xaxis);
  }
  kt1 = kt1 * (rkt * cos(phi) / kt1.P());

  Vec4 kt2cms = Q.Boost(pijt);
  kt2cms = kt2cms.Cross(kt1);
  kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.P());
  Vec4 kt2 = Q.BoostBack(kt2cms);

  // Conversion to {i, j, k} basis
  Vec4 pi = pijt * z + pkt * ((1. - z) * y) + kt1 + kt2;
  Vec4 pj = pijt * (1. - z) + pkt * (z * y) - kt1 - kt2;
  Vec4 pk = pkt * (1. - y);

  // No need to do *kinematics[0], for arrays the elements are already
  // pointers
  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
}

#endif  // SHOWER_KINEMATICS_CUH_