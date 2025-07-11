#include "shower.h"

void shower::make_kinematics(vec4* kinematics, const double z, const double y,
                             const double phi, const vec4 pijt, const vec4 pkt,
                             int sf) const {
  /**
   * @brief Generate the post-emission kinematics for the given splitting
   *
   * @param kinematics: Array to store the post-emission momenta
   * @param z: Splitting variable
   * @param y: y(t,z) value
   * @param phi: Azimuthal angle
   * @param pijt: Emitter momentum
   * @param pkt: Spectator momentum
   * @param sf: Splitting function code
   */

  // --------------------------------------------------
  // k_perp Calculation - Rodrigues' Rotation Formula

  // Splitter and Spectator
  vec4 m = pijt + pkt;

  // generating the transverse momentum
  double pt = sqrt(z * (1. - z) * y * m.m2());

  // Kt vector in the dipole frame
  vec4 kt_dipole = vec4(0., pt * cos(phi), pt * sin(phi), 0.);

  // Rotate to frame where pijt is along the z-axis
  vec4 z_axis(0., 0., 0., 1.);
  vec4 pijt_b2b = m.boost(pijt);
  pijt_b2b = pijt_b2b * (1. / pijt_b2b.p());
  double theta = acos(z_axis.dot(pijt_b2b));
  vec4 axis = z_axis.cross(pijt_b2b);
  if (axis.p() > 1e-6) {
    axis = axis * (1. / axis.p());
  }
  vec4 kt_b2b = kt_dipole;
  // No need to rotate if the angle is zero
  if (!((theta < 1e-6) || (axis.p() < 1e-6))) {
    // Rodrigues' Rotation Formula
    // v = vcosT + (k x v)sinT + k(k.v)(1 - cosT)
    kt_b2b = kt_dipole * cos(theta) + axis.cross(kt_dipole) * sin(theta) +
             axis * (axis.dot(kt_dipole)) * (1. - cos(theta));
  }

  // Boost back to the original frame
  vec4 kt = m.boost_back(kt_b2b);

  // --------------------------------------------------
  // k_perp Calculation - Gram-Schmidt Process (OLD)

  // vec4 q = pijt + pkt;
  // double rkt = sqrt(calculate_t(z, y, q.m2(), sf));

  // vec4 kt1 = pijt.cross(pkt);
  // if (kt1.p() < 1.e-6) {
  //   vec4 xaxis(0., 1., 0., 0.);
  //   kt1 = pijt.cross(xaxis);
  // }
  // kt1 = kt1 * (rkt * cos(phi) / kt1.p());

  // vec4 kt2cms = q.boost(pijt);
  // kt2cms = kt2cms.cross(kt1);
  // kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.p());
  // vec4 kt2 = q.boost_back(kt2cms);

  // vec4 kt = kt1 + kt2;

  // --------------------------------------------------
  // Post emission momenta by splitting type

  vec4 pi, pj, pk;

  pi = pijt * z + pkt * ((1. - z) * y) + kt;
  pj = pijt * (1. - z) + pkt * (z * y) - kt;
  pk = pkt * (1. - y);

  // --------------------------------------------------
  // Set the kinematics

  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
  kinematics[3] = pijt;
  kinematics[4] = pkt;
  kinematics[5] = kt;
}