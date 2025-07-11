#include "matrix.h"

void matrix::lep_lo(event& ev) {
  /**
   * @brief generate the leading order interaction e+ e- -> q qbar
   *
   * The Leading order differential cross section for e+ e- -> q qbar is
   * given by:
   *
   * dsigma = drho1 drho2 * 1/(2 s_hat) * 1/(8 pi) * N_C * 1/4 Sum_spins |M|^2
   *
   * @param ev the event object
   */

  // Determine the flavour of the quark and antiquark
  int fl = static_cast<int>(ev.gen_random() * 5) + 1;

  // Two random numbers for cos(theta) and phi
  double rho_1 = ev.gen_random();
  double rho_2 = ev.gen_random();

  // Generate the Phase Space Variables cos(theta), phi
  double ct = -1. + 2. * rho_1;
  double phi = 2. * M_PI * rho_2;

  // Kinematics
  double ecms = sqrt(s_hat);
  double p0 = 0.5 * ecms;               // 0.5 * Ecms
  double st = std::sqrt(1. - ct * ct);  // sin(theta)
  vec4 pa = vec4(p0, 0., 0., p0);       // e+
  vec4 pb = vec4(p0, 0., 0., -p0);      // e-
  vec4 p1 = vec4(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);     // q
  vec4 p2 = vec4(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);  // qb

  // Generate the particles
  particle p[4] = {particle(), particle(), particle(), particle()};
  p[0] = particle(-11, pa, 0, 0, 1000.);  // dummy eta
  p[1] = particle(11, pb, 0, 0, 1000.);   // dummy eta
  p[2] = particle(fl, p1, 1, 0);
  p[3] = particle(-fl, p2, 0, 1);

  // Calculate the matrix element
  double lome = me2(fl, (pa + pb).m2(), (pa - p1).m2());
  lome *= k_nc;  // Three possible colour states

  // Calculate the cross-section
  double dxs = (1. / (2. * s_hat)) * (1. / (8. * M_PI)) * lome;
  dxs *= static_cast<double>(k_nf) * GeV_minus_2_to_pb;  // 5 flavours + units

  // Set the particles
  for (int i = 0; i < 4; i++) {
    ev.set_particle(i, p[i]);
  }
  ev.set_hard(4);

  // Store the Matrix Element and Cross Section
  ev.set_me2(lome);
  ev.set_dxs(dxs);
}