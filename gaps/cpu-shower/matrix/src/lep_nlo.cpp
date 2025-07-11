#include "matrix.h"

void matrix::lep_nlo(event &ev) {
  /**
   * @brief generate the next-to-leading order interaction e+ e- -> q qbar g
   *
   * There are two types of events possible
   * - H Event: There is an emission in the LO system
   * - S Evennt: The LO and Virtual Contributions
   *
   * @params ev the event object
   */

  // To prevent duplicate code, generate LO event first
  // we'll use the particles and the me2 from this
  lep_lo(ev);

  // define common variables
  double me2_lo = ev.get_me2();
  double dxs_lo = ev.get_dxs();
  double dxs_nlo;

  // Do ee -> qqbar g based on the probability ws (1/4 events)
  if (ev.gen_random() < ws) {
    /**
     * Catani Seymour Real Corrections
     * -------------------------------
     *
     * Like the dipole shower, here we add a new parton to the event. We do this
     * using the catani seymour dipole kinematics. The dxs is also adjusted by
     * a factor. See Catani-Seymour Appendix D.7 or Black Book pages 150 -> 154.
     *
     * The H event cross section is given by:
     * dxs = drho1 drho2 drho3 drho4 drho5 * 1/(2 s_hat) * 1/(8 pi)
     *       * s_hat / (16 pi^2) * (|M|^2 - D132 - D231) * (1 - y)
     */

    // -------------------------------------------------------------------------
    // Generate the Emission

    // Need three random numbers for y, z and phi
    double rho_3 = ev.gen_random();
    double rho_4 = ev.gen_random();
    double rho_5 = ev.gen_random();

    // Generate y, z and phi for the emission
    double y = rho_3;  // std::pow(ev.gen_random(), 1. / (1. - ye));
    double z = rho_4;  // std::pow(ev.gen_random(), 1. / (1. - ze));
    double phi = 2. * M_PI * rho_5;

    // Randomly determine which parton is the emitter
    vec4 pijt, pkt;
    bool ij_is_quark = ev.gen_random() < 0.5;
    if (ij_is_quark) {
      pijt = ev.get_particle(2).get_mom();
      pkt = ev.get_particle(3).get_mom();
    } else {
      pijt = ev.get_particle(3).get_mom();
      pkt = ev.get_particle(2).get_mom();
    }

    // The SAME Kinematics as the dipole shower
    vec4 m = pijt + pkt;
    double rkt = sqrt(m.m2() * y * z * (1. - z));
    vec4 kt1 = pijt.cross(pkt);
    if (kt1.p() < 1.e-6) {
      vec4 xaxis(0., 1., 0., 0.);
      kt1 = pijt.cross(xaxis);
    }
    kt1 = kt1 * (rkt * cos(phi) / kt1.p());
    vec4 kt2cms = m.boost(pijt);
    kt2cms = kt2cms.cross(kt1);
    kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.p());
    vec4 kt2 = m.boost_back(kt2cms);
    vec4 pi = pijt * z + pkt * ((1. - z) * y) + kt1 + kt2;
    vec4 pj = pijt * (1. - z) + pkt * (z * y) - kt1 - kt2;
    vec4 pk = pkt * (1. - y);

    // -------------------------------------------------------------------------
    // Adjust the Cross Section

    // Define Vectors for quark(1), antiquark(2), gluon(3)
    vec4 p1 = ij_is_quark ? pi : pk;
    vec4 p2 = ij_is_quark ? pk : pi;
    vec4 p3 = pj;
    double s12 = (p1 + p2).m2();
    double s13 = (p1 + p3).m2();
    double s23 = (p2 + p3).m2();
    double Q2 = s12 + s13 + s23;
    int fl = ev.get_particle(2).get_pid();

    // Real Emission Term
    double real =
        me2_lo / Q2 * (s23 / s13 + s13 / s23 + 2. * s12 * Q2 / (s13 * s23));
    real *= k_cf * 8 * M_PI * as(s_hat);

    // Subtraction term - quark emitter
    double z1 = s12 / (s12 + s23);
    double y132 = s13 / (s12 + s13 + s23);
    double D132 = 1. / s13 * (2. / (1. - z1 * (1. - y132)) - (1. + z1));
    vec4 p13t = p1 + p3 - p2 * (y132 / (1. - y132));
    D132 *= me2(fl, Q2, (ev.get_particle(0).get_mom() - p13t).m2()) * k_nc;
    D132 *= k_cf * 8 * M_PI * as(s_hat);

    // Subtraction term - antiquark emitter
    double z2 = s12 / (s12 + s13);
    double y231 = s23 / (s12 + s13 + s23);
    double D231 = 1. / s23 * (2. / (1. - z2 * (1. - y231)) - (1. + z2));
    vec4 p23t = p2 + p3 - p1 * (y231 / (1. - y231));
    D231 *= me2(fl, Q2, (ev.get_particle(1).get_mom() - p23t).m2()) * k_nc;
    D231 *= k_cf * 8 * M_PI * as(s_hat);

    // Veto very small virtualities - no subtraction
    if (y132 < amin || y231 < amin) {
      ev.set_dxs(0.);
      return;
    }

    // Calculate the cross-section
    dxs_nlo = (1. / (2. * s_hat)) * (1. / (8. * M_PI));
    dxs_nlo *= s_hat / (16. * M_PI * M_PI);
    dxs_nlo *= (real - D132 - D231);
    dxs_nlo *= (1. - y);
    dxs_nlo *= static_cast<double>(k_nf) * GeV_minus_2_to_pb;

    // Adjust dxs for the nlo weighting
    dxs_nlo /= ws;

    // -------------------------------------------------------------------------
    // Finalise the Event

    // Create the gluon particle
    particle g;
    g.set_pid(21);

    // set the momenta to the relevant partons
    if (ij_is_quark) {
      // quark
      ev.set_particle_mom(2, pi);
      ev.set_particle_col(2, 2);
      // gluon
      g.set_mom(pj);
      g.set_col(1);
      g.set_acol(2);
      // antiquark
      ev.set_particle_mom(3, pk);
    } else {
      // antiquark
      ev.set_particle_mom(3, pi);
      ev.set_particle_acol(3, 2);
      // gluon
      g.set_mom(pj);
      g.set_col(2);
      g.set_acol(1);
      // quark
      ev.set_particle_mom(2, pk);
    }

    // Add the gluon to the event
    ev.set_particle(4, g);
    ev.set_hard(5);

    // Set the new cross-section
    ev.set_dxs(dxs_nlo);

    // -------------------------------------------------------------------------
    // Matching - Set the starting scale for the shower

    // This is t of the first emission
    double t = ij_is_quark ? s13 * z1 * (1. - z1) : s23 * z2 * (1. - z2);

    ev.set_shower_t(t);
    ev.set_shower_c(2);

    return;
  }

  // Do ee -> qqbar (Born + Virtual) for the rest of the events
  else {
    /**
     * Catani Seymour Virtual Corrections
     * -----------------------------------
     *
     * For this process, the virtual correction can be calculated analytically,
     * so here we just make a note, but just put in the end result.
     *
     * We will not be using the extra term L = mu^2/s as we set mu^2 = s:
     * L = m.log(mu*mu/(lo[0][2].mom+lo[0][3].mom).M2())
     *
     * The Divergent terms are:
     * V = - 2/e^2 - 3/e - 8 + pi^2
     *
     * We Add two Catani-Seymour counterterms:
     * I(132) = 1/e^2 + 3/2e + 5 - pi^2/2
     * I(231) = 1/e^2 + 3/2e + 5 - pi^2/2
     *
     * The final result is:
     * V + I(132) + I(231) = 2
     * Coeff = CF as(mu^2) / 2pi
     * Factor = Coeff * (V + I(132) + I(231)) = CF as(mu^2) / pi
     *
     * You mutliply the LO dsigma by (1 + factor) to get the NLO dsigma
     */

    // -------------------------------------------------------------------------
    // Calculate the NLO cross-section, params[0] = ecms

    // Calculate the factor to adjust the cross-section
    double factor = k_cf * as(s_hat) / M_PI;
    dxs_nlo = dxs_lo * (1. + factor);

    // Adjust dxs for the nlo weighting ws
    dxs_nlo /= (1. - ws);

    // Set the new cross-section
    ev.set_dxs(dxs_nlo);

    // -------------------------------------------------------------------------
    // Matching - Set the starting scale for the shower

    double sij =
        (ev.get_particle(2).get_mom() + ev.get_particle(3).get_mom()).m2();

    ev.set_shower_t(sij);
    ev.set_shower_c(1);

    return;
  }
}