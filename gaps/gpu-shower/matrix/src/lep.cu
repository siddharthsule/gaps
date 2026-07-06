#include "matrix.cuh"

__device__ double matrix::me2_ee2Zy2qq(int fl, double s, double t) const {
  /**
   * @brief Generate the matrix element squared for massless 2 -> 2 scattering
   * of e+ e- -> q qbar
   *
   * matrix element squared for massess 2x2 scattering via a virtual photon
   * or a z boson. see peskin and schroeder page 156 and pink book page 54. This
   * contains both the photon and the z boson contributions. The matrix element
   * squared is given by:
   *
   * me2 = 16 pi^2 alpha^2 ( (1+cos^2theta)(3qf + ...) + costheta(4qf...) )
   *
   * @param fl The flavour of the quark
   * @param s The Mandelstam s
   * @param t The Mandelstam t
   * @return double: The matrix element squared
   */

  // constants Ae, Ve, Af, Vf, Qe, Qf
  double qe = -1.;
  double ae = -0.5;
  double ve = ae - 2. * qe * sin2tw;
  double qf = (abs(fl) == 2 || abs(fl) == 4) ? 2. / 3. : -1. / 3.;
  double af = (abs(fl) == 2 || abs(fl) == 4) ? 0.5 : -0.5;
  double vf = af - 2. * qf * sin2tw;

  // Kappa and Chi Functions (Breit-Wigner Here!)
  double kappa = 1. / (4. * sin2tw * (1. - sin2tw));
  double chi1 = kappa * s * (s - mz2) / (pow(s - mz2, 2.) + gz2 * mz2);
  double chi2 = pow(kappa * s, 2.) / (pow(s - mz2, 2.) + gz2 * mz2);

  // Scattering Angle of CoM Final State Fermions
  // Can write in terms of s and t, but this is easier to read
  double cos_theta = 1. + 2. * t / s;

  // Slight difference from Pink Book pg 54: book assumes qe = -1 so skips
  double term1 = (1. + pow(cos_theta, 2.)) *
                 (pow(qf * qe, 2.) + 2. * (qf * qe * vf * ve) * chi1 +
                  (ae * ae + ve * ve) * (af * af + vf * vf) * chi2);

  double term2 = cos_theta * (4. * qe * qf * ae * af * chi1 +
                              8. * ae * ve * af * vf * chi2);

  // Output: 16 pi^2 aEM^2 ( (1+cos^2theta)(3qf + ...) + costheta(4qf...) )
  return pow(4. * M_PI * alpha, 2.) * (term1 + term2);
}

__global__ void lep_lo_kernel(matrix* matrix, event* events, int n) {
  /**
   * @brief generate the leading order interaction e+ e- -> q qbar
   *
   * The Leading order differential cross section for e+ e- -> q qbar is
   * given by:
   *
   * dsigma = drho1 drho2 * 1/(2s) * 1/(8 pi) * N_C * 1/4 Sum_spins |M|^2
   *
   * @param matrix matrix element generator
   * @param events array of event records
   * @param n number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  // ---------------------------------------------

  // Determine the flavour of the quark and antiquark
  int fl = static_cast<int>(ev.gen_random() * 5) + 1;

  // Two random numbers for cos(theta) and phi
  double rho_1 = ev.gen_random();
  double rho_2 = ev.gen_random();

  // Generate the Phase Space Variables cos(theta), phi
  double ct = -1. + 2. * rho_1;
  double phi = 2. * M_PI * rho_2;

  // Kinematics
  double ecms = matrix->root_s;     // Ecms
  double p0 = 0.5 * ecms;           // 0.5 * Ecms
  double st = sqrt(1. - ct * ct);   // sin(theta)
  vec4 pa = vec4(p0, 0., 0., p0);   // e+
  vec4 pb = vec4(p0, 0., 0., -p0);  // e-
  vec4 p1 = vec4(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);     // q
  vec4 p2 = vec4(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);  // qb

  // Generate the particles
  particle p[4] = {particle(), particle(), particle(), particle()};
  p[0] = particle(-11, pa, 0, 0, 1000.);  // dummy eta
  p[1] = particle(11, pb, 0, 0, 1000.);   // dummy eta
  p[2] = particle(fl, p1, 1, 0);
  p[3] = particle(-fl, p2, 0, 1);

  // Calculate the matrix element
  double lome = matrix->me2_ee2Zy2qq(fl, (pa + pb).m2(), (pa - p1).m2());
  lome *= k_nc;  // Three possible colour states

  // Calculate the cross-section
  double dxs =
      (1. / (2. * matrix->root_s * matrix->root_s)) * (1. / (8. * M_PI)) * lome;
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

void lep_lo(thrust::device_vector<event>& d_events, matrix* matrix, int blocks,
            int threads) {
  /**
   * @brief Wrapper for the leading order LEP kernel
   *
   * @param d_events device vector of events
   * @param matrix matrix element generator
   * @param blocks number of blocks
   * @param threads number of threads
   */

  debug_msg("running @lep_lo_kernel");
  lep_lo_kernel<<<blocks, threads>>>(
      matrix, thrust::raw_pointer_cast(d_events.data()), d_events.size());
  sync_gpu_and_check("lep_lo_kernel");
}

__global__ void lep_nlo_kernel(matrix* matrix, alpha_s* as, event* events,
                               int n) {
  /**
   * @brief generate the next-to-leading order interaction e+ e- -> q qbar g
   *
   * There are two types of events possible
   * - H Event: There is an emission in the LO system
   * - S Event: The LO and Virtual Contributions
   *
   * @param matrix matrix element generator
   * @param as alpha_s generator
   * @param events array of event records
   * @param n number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  // ---------------------------------------------

  // define common variables
  double me2_lo = ev.get_me2();
  double dxs_lo = ev.get_dxs();
  double dxs_nlo;

  // Do ee -> qqbar g based on the probability ws (1/4 events)
  if (ev.gen_random() < matrix->ws) {
    /**
     * Catani Seymour Real Corrections
     * -------------------------------
     *
     * Like the dipole shower, here we add a new parton to the event. We do this
     * using the catani seymour dipole kinematics. The dxs is also adjusted by
     * a factor. See Catani-Seymour Appendix D.7 or Black Book pages 150 -> 154.
     *
     * The H event cross section is given by:
     * dsigma = drho1 drho2 drho3 drho4 drho5 * 1/(2 s) * 1/(8 pi)
     *          * s / (16 pi^2) * (|M|^2 - D132 - D231) * (1 - y)
     */

    // -------------------------------------------------------------------------
    // Generate the Emission

    // Need three random numbers for y, z and phi
    double rho_3 = ev.gen_random();
    double rho_4 = ev.gen_random();
    double rho_5 = ev.gen_random();

    // Generate y, z and phi for the emission
    double y = rho_3;
    double z = rho_4;
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
    real *= k_cf * 8 * M_PI * (*as)(matrix->root_s * matrix->root_s);

    // Subtraction term - quark emitter
    double z1 = s12 / (s12 + s23);
    double y132 = s13 / (s12 + s13 + s23);
    double D132 = 1. / s13 * (2. / (1. - z1 * (1. - y132)) - (1. + z1));
    vec4 p13t = p1 + p3 - p2 * (y132 / (1. - y132));
    D132 *= matrix->me2_ee2Zy2qq(fl, Q2,
                                 (ev.get_particle(0).get_mom() - p13t).m2()) *
            k_nc;
    D132 *= k_cf * 8 * M_PI * (*as)(matrix->root_s * matrix->root_s);

    // Subtraction term - antiquark emitter
    double z2 = s12 / (s12 + s13);
    double y231 = s23 / (s12 + s13 + s23);
    double D231 = 1. / s23 * (2. / (1. - z2 * (1. - y231)) - (1. + z2));
    vec4 p23t = p2 + p3 - p1 * (y231 / (1. - y231));
    D231 *= matrix->me2_ee2Zy2qq(fl, Q2,
                                 (ev.get_particle(1).get_mom() - p23t).m2()) *
            k_nc;
    D231 *= k_cf * 8 * M_PI * (*as)(matrix->root_s * matrix->root_s);

    // Veto very small virtualities - no subtraction
    double ymin = 1e-10;
    if (y132 < ymin || y231 < ymin) {
      ev.set_dxs(0.);
      return;
    }

    // Calculate the cross-section
    dxs_nlo =
        (1. / (2. * matrix->root_s * matrix->root_s)) * (1. / (8. * M_PI));
    dxs_nlo *= matrix->root_s * matrix->root_s / (16. * M_PI * M_PI);
    dxs_nlo *= (real - D132 - D231);
    dxs_nlo *= (1. - y);
    dxs_nlo *= static_cast<double>(k_nf) * GeV_minus_2_to_pb;

    // Adjust dxs for the nlo weighting
    dxs_nlo /= matrix->ws;

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
    ev.add_emission(g);

    // Set the new cross-section
    ev.set_dxs(dxs_nlo);

    // -------------------------------------------------------------------------
    // Matching - Set the starting scale for the shower

    // This is t of the first emission
    double t =
        ij_is_quark ? s13 * y132 * z1 * (1. - z1) : s23 * y231 * z2 * (1. - z2);

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
    double factor = k_cf * (*as)(matrix->root_s * matrix->root_s) / M_PI;
    dxs_nlo = dxs_lo * (1. + factor);

    // Adjust dxs for the nlo weighting ws
    dxs_nlo /= (1. - matrix->ws);

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

void lep_nlo(thrust::device_vector<event>& d_events, matrix* matrix,
             alpha_s* as, int blocks, int threads) {
  /**
   * @brief wrapper for the lep_nlo kernel
   *
   * @param d_events device vector of events
   * @param matrix matrix element generator
   * @param as alpha_s generator
   * @param blocks number of blocks
   * @param threads number of threads
   */
  debug_msg("running @lep_nlo_kernel");
  lep_nlo_kernel<<<blocks, threads>>>(
      matrix, as, thrust::raw_pointer_cast(d_events.data()), d_events.size());
  sync_gpu_and_check("lep_nlo_kernel");
}