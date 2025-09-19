#include "shower.h"

double shower::calculate_y(double t, double z, double sijk, int sf) const {
  /**
   * @brief Calculate y(z,t) for the given splitting
   *
   * @param t: Evolution variable
   * @param z: Splitting variable
   * @param sijk: Invariant mass of the emitter and spectator
   * @param sf: Splitting function code
   * @return double: y(t,z)
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  switch (splitting_case) {
    // FF splittings: z and y
    case 0:
      return t / sijk / z / (1. - z);
      break;

    // FI splittings: z and x (here, y)
    case 1:
      return 1 / ((t / sijk / z / (1. - z)) + 1);
      break;

    // IF splittings: x (here, z) and u (here, y)
    case 2:
      return 0.5 * (1 - sqrt(1 - (4 * t * z) / (sijk * (1 - z))));
      break;

    // II splittings: x (here, z) and v (here, y)
    case 3:
      return (1 - z) / 2. *
             (1 - sqrt(1 - (4 * t * z) / (sijk * (1 - z) * (1 - z))));
      break;
  }
  return 0.;
}

double shower::calculate_t(double z, double y, double sijk, int sf) const {
  /**
   * @brief Calculate t(z,y) for the given splitting
   *
   * @param z: Splitting variable
   * @param y: Momentum fraction of the spectator
   * @param sijk: Invariant mass of the emitter and spectator
   * @param sf: Splitting function code
   * @return double: t value
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  switch (splitting_case) {
    // FF splittings: z and y
    case 0:
      return z * (1. - z) * y * sijk;
      break;

    // FI splittings: z and x (here, y)
    case 1:
      return z * (1 - z) * (1 - y) / y * sijk;
      break;

    // IF splittings: x (here, z) and u (here, y)
    case 2:
      return y * (1 - y) * (1 - z) / z * sijk;
      break;

    // II splittings: x (here, z) and v (here, y)
    case 3:
      return (1 - z - y) * y / z * sijk;
      break;
  }
  return 0.;
}

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
  double pt = sqrt(calculate_t(z, y, m.m2(), sf));

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

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  switch (splitting_case) {
    // FF splittings: z and y
    case 0:
      pi = pijt * z + pkt * ((1. - z) * y) + kt;
      pj = pijt * (1. - z) + pkt * (z * y) - kt;
      pk = pkt * (1. - y);
      break;

    // FI splittings: z and x (here, y)
    case 1:
      pi = pijt * z + pkt * ((1. - z) * (1. - y) / y) + kt;
      pj = pijt * (1. - z) + pkt * (z * (1. - y) / y) - kt;
      pk = pkt * (1. / y);
      break;

    // IF splittings: x (here, z) and u (here, y) - Boost after emission
    case 2:
      pi = pijt * (1. / z);
      pj = pijt * ((1. - y) * (1. - z) / z) + pkt * y - kt;
      pk = pijt * (y * (1. - z) / z) + pkt * (1. - y) + kt;
      /*
      // Alternative Map - NOT USING THIS
      pi = pijt * ((1 - y) / (z - y)) + pkt * ((y / z) * (1. - z) / (z - y)) +
           kt * (1. / (y - z));
      pj = pijt * ((1 - z) / (z - y)) + pkt * ((y / z) * (1. - y) / (z - y)) +
           kt * (1. / (y - z));
      pk = pkt * ((z - y) / z);
      */
      break;

    // II splittings: x (here, z) and v (here, y) - Boost after emission
    case 3:
      pi = pijt * (1. / z);
      pj = pijt * ((1. - z - y) / (z)) + pkt * y - kt;
      pk = pkt;
      break;
  }

  // --------------------------------------------------
  // Set the kinematics

  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
  kinematics[3] = pijt;  // for IF and II Boost after emission
  kinematics[4] = pkt;
  kinematics[5] = kt;
}

// -----------------------------------------------------------------------------
// Boost after emission

void shower::ii_boost_after_emission(event& ev, vec4* kinematics) const {
  /**
   * @brief Boost the momenta after the emission for II splittings
   *
   * @param ev: Event object
   * @param kinematics: Array of post-emission momenta
   */

  vec4 ktil = kinematics[3] + kinematics[4];  // pij + pk OR  pijt + pkt
  vec4 k = kinematics[0] - kinematics[1] + kinematics[2];  // pi - pj + pk

  // temporary momenta for the boost
  for (int i = 0; i < ev.get_size(); i++) {
    // Avoid initial state particles
    if (ev.get_particle(i).is_initial()) {
      continue;
    }

    // Avoid the emission - should be the last element
    if (i == ev.get_size() - 1) {
      continue;
    }

    // Boost the original momenta
    vec4 p = ev.get_particle(i).get_mom();

    vec4 q1 = p;
    vec4 q2 = (ktil + k) * (-2.) * (1. / (ktil + k).m2()) * ((ktil + k) * p);
    vec4 q3 = k * (2.) * (1. / k.m2()) * (ktil * p);
    vec4 q = q1 + q2 + q3;

    // Assign the new momenta
    ev.set_particle_mom(i, q);
  }
}