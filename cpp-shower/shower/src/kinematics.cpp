#include "shower.h"

void shower::make_kinematics(vec4* kinematics, const double z, const double y,
                             const double phi, const vec4 pijt,
                             const vec4 pkt) {
  vec4 q = pijt + pkt;

  // generating the momentum (0, kt1, kt2, 0)
  double rkt = sqrt(q.m2() * y * z * (1. - z));

  vec4 kt1 = pijt.cross(pkt);
  if (kt1.p() < 1.e-6) {
    vec4 xaxis(0., 1., 0., 0.);
    kt1 = pijt.cross(xaxis);
  }
  kt1 = kt1 * (rkt * cos(phi) / kt1.p());

  vec4 kt2cms = q.boost(pijt);
  kt2cms = kt2cms.cross(kt1);
  kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.p());
  vec4 kt2 = q.boost_back(kt2cms);

  // conversion to {i, j, k} basis
  vec4 pi = pijt * z + pkt * ((1. - z) * y) + kt1 + kt2;
  vec4 pj = pijt * (1. - z) + pkt * (z * y) - kt1 - kt2;
  vec4 pk = pkt * (1. - y);

  // no need to do *kinematics[0], for arrays the elements are already pointers
  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
}