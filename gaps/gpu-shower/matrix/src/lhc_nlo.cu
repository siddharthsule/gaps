#include "matrix.cuh"

// -----------------------------------------------------------------------------
// Emulating DiLog

__device__ double dilogarithm(double x) {
  /**
   * @brief calculate the dilogarithm of x
   *
   * Li2(x) = Sum_{n=1}^{inf} (x^n)/(n^2)
   *
   * We do this for until a precision of 1e-10 is reached
   *
   * @param x the value to calculate the dilogarithm of
   * @return double the dilogarithm of x
   */

  double sum = 0;
  for (int n = 1; n < 100; n++) {
    double temp = pow(x, n) / (n * n);
    if (fabs(temp) < 1e-12) {
      break;
    }
    sum += temp;
  }

  return sum;
}

// -----------------------------------------------------------------------------
// Matrix Element Calculations

__device__ double matrix::me2qqZ(int fl, double s) const {
  /**
   * @brief Generate the matrix element squared for just q qbar -> Z
   *
   * @param fl the flavour of the quark
   * @param s the Mandelstam s variable
   */

  // For Comparison with Herwig
  // alpha = 1. / 128.91;
  // sin2tw = 0.232;

  // constants A, V, Q
  double q = (abs(fl) == 2 || abs(fl) == 4) ? 2. / 3. : -1. / 3.;
  double a = (abs(fl) == 2 || abs(fl) == 4) ? 0.5 : -0.5;
  double v = a - 2. * q * sin2tw;

  // Fermi Constant \sqrt{2} G_F
  double root2_gf = (4. * M_PI * alpha) / (mz * mz * sin2tw * (1. - sin2tw));

  // Calculate the Matrix Element
  double me2;
  me2 = root2_gf;            // Root 2 G_F
  me2 *= mz * mz * mz * mz;  // M_Z^4
  me2 *= (a * a + v * v);    // Couplings
  me2 *= (1. / 4.);          // spin average
  me2 *= 1 / k_nc;           // Three possible initial colour states

  return me2;
}

// -----------------------------------------------------------------------------
// Utility Kernels

__global__ void multiply_dxs_by_pdf(event *events, int n, double *xf) {
  /**
   * @brief Multiply the cross section of each event by the PDFs
   *
   * @param events array of event records
   * @param n number of events
   * @param xf_a array of PDF values
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event &ev = events[idx];
  // ---------------------------------------------

  // Mutliply by PDF
  double pdf = xf[idx];
  ev.set_dxs(ev.get_dxs() * pdf);
}

// -----------------------------------------------------------------------------
// LO Event Generation - pp -> Z

__global__ void lo_event(event *events, int n, matrix *matrix, int *fl_a,
                         int *fl_b, double *xa, double *xb, double *q2) {
  /**
   * @brief LO event generation for pp -> Z
   *
   * @param events array of event records
   * @param n number of events
   * @param matrix pointer to the matrix object
   * @param fl_a array of flavour indices for parton A
   * @param fl_b array of flavour indices for parton B
   * @param xa array of momentum fractions for parton A
   * @param xb array of momentum fractions for parton B
   * @param q2 array of invariant mass squared for the events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event &ev = events[idx];
  // ---------------------------------------------

  // We generate an On Shell Z boson, so we cannot run the lo function and use
  // output, we must generate the event ourselves (and fix s_hat)

  // LHC: Bias Flavout based on contribution (based on cross section)
  // Flavours:    [down,       up,         strange,    charm,      bottom    ]
  // Percentages: [0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279]
  // Cumulatives: [0.36677679, 0.80186718, 0.91861688, 0.97129721, 1.00000000]
  int fl;
  double pd[5] = {0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279};
  double pc[5] = {0.36677679, 0.80186718, 0.91861688, 0.97129721, 1.00000000};
  // double pd[5] = {1., 0., 0., 0., 0.};
  // double pc[5] = {1., 1., 1., 1., 1.};
  double r = ev.gen_random();
  for (int i = 0; i < 5; ++i) {
    if (r < pc[i]) {
      fl = i + 1;  // Flavors are 1-based
      break;
    }
  }

  // On Shell Boson, so s_hat = mz2
  double s_hat = mz * mz;

  // So only generate one random number!
  double rho_1 = ev.gen_random();

  // generate y
  double lim = fmin(100., 0.5 * log((matrix->root_s * matrix->root_s) / s_hat));
  double y = -lim + 2. * lim * rho_1;

  // Momentum Fractions
  double eta_a = sqrt(s_hat) / matrix->root_s * exp(y);
  double eta_b = sqrt(s_hat) / matrix->root_s * exp(-y);

  // Quark / AntiQuark Momentum
  double p0 = 0.5 * matrix->root_s;
  vec4 pa = vec4(eta_a * p0, 0., 0., eta_a * p0);
  vec4 pb = vec4(eta_b * p0, 0., 0., -eta_b * p0);

  // Two Possibilities: P(q) P(qbar) <-> P(qbar) P(q)
  // In this case, swap momentum and momentum fraction
  if (ev.gen_random() < 0.5) {
    vec4 temp_p = pa;
    pa = pb;
    pb = temp_p;
    double temp_eta = eta_a;
    eta_a = eta_b;
    eta_b = temp_eta;
  }

  // Z Momentum
  vec4 pz = pa + pb;

  // Generate the particles
  particle p[3] = {particle(), particle(), particle()};
  p[0] = particle(fl, pa, 0, 1, eta_a);
  p[1] = particle(-fl, pb, 1, 0, eta_b);
  p[2] = particle(23, pz, 0, 0);

  // Calculate the Matrix Element Squared
  double lome = matrix->me2qqZ(fl, (pa + pb).m2());

  // Store PDF Calculation details
  fl_a[idx] = fl;
  fl_b[idx] = -fl;
  xa[idx] = eta_a;
  xb[idx] = eta_b;
  q2[idx] = s_hat;

  // Calculate the Cross Section
  double dxs;
  dxs = 2.;  // Two Possible Orientations (P(q) P(qbar) or P(qbar) P(q))
  dxs *= M_PI;
  dxs *= log((matrix->root_s * matrix->root_s) / s_hat);  // y Sampling
  dxs *= lome;                                            // ME2
  dxs *= (1. / pow(s_hat, 2));                            // Leftover
  dxs *= GeV_minus_2_to_pb;                               // units
  dxs /= pd[abs(fl) - 1];                                 // Flavour Selection

  // Set the particles
  for (int i = 0; i < 3; i++) {
    ev.set_particle(i, p[i]);
  }
  ev.set_hard(3);

  // Store the Matrix Element and Cross Section
  ev.set_me2(lome);
  ev.set_dxs(dxs);

  // LO Shower Settings
  ev.set_shower_t(pz.m2());
  ev.set_shower_c(1);

  return;
}

// -----------------------------------------------------------------------------
// H Event Generation - pp -> Zj

__global__ void h_event(event *events, int n, matrix *matrix, alpha_s *as,
                        int *fl_a, int *fl_b, double *xa, double *xb,
                        double *q2pdf) {
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event &ev = events[idx];
  // ---------------------------------------------

  // ---------------------------------------------------------------------------
  // Is this an H Event? If not, return
  if (ev.gen_random() > matrix->ws) {
    return;
  }

  // ---------------------------------------------------------------------------
  // Need to do this again because new kernel

  // Get the matrix element squared
  double me2_lo = ev.get_me2();
  double dxs_lo = ev.get_dxs();

  // Momentum Fractions
  double eta_q = ev.get_particle(0).get_eta();
  double eta_qbar = ev.get_particle(1).get_eta();

  // Scale Choice = mz2, the only final state particle in Born!
  vec4 pz = ev.get_particle(2).get_mom();
  double mu2 = pz.m2();

  // flavour of the qqbar pair
  int fl = abs(ev.get_particle(0).get_pid());

  // ---------------------------------------------------------------------------
  // Now onto the H Event Generation

  // Setup: Choose the emitter
  bool ij_is_quark = ev.gen_random() < 0.5;
  bool is_q2qg = ev.gen_random() < 0.5;

  // Need three more random numbers
  double rho_2 = ev.gen_random();
  double rho_3 = ev.gen_random();
  double rho_4 = ev.gen_random();

  // Generate x, v and phi for the emission
  double x_min = ij_is_quark ? eta_q : eta_qbar;
  double x = x_min + (1. - x_min) * rho_2;
  double v = (1. - x) * rho_3;
  double phi = 2. * M_PI * rho_4;

  // Kinematics
  vec4 pait, pbt;
  if (ij_is_quark) {
    pait = ev.get_particle(0).get_mom();
    pbt = ev.get_particle(1).get_mom();
  } else {
    pait = ev.get_particle(1).get_mom();
    pbt = ev.get_particle(0).get_mom();
  }

  // The SAME Kinematics as the dipole shower
  vec4 m = pait + pbt;
  double rkt = sqrt((1 - x - v) * v / x * m.m2());
  vec4 kt1 = pait.cross(pbt);
  if (kt1.p() < 1.e-6) {
    vec4 xaxis(0., 1., 0., 0.);
    kt1 = pait.cross(xaxis);
  }
  kt1 = kt1 * (rkt * cos(phi) / kt1.p());
  vec4 kt2cms = m.boost(pait);
  kt2cms = kt2cms.cross(kt1);
  kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.p());
  vec4 kt2 = m.boost_back(kt2cms);
  vec4 kt = kt1 + kt2;

  // Calculate the momenta
  vec4 pa = pait * (1. / x);
  vec4 pi = pait * ((1. - x - v) / (x)) + pbt * v + kt;
  vec4 pb = pbt;

  // Boost the Z Boson Momenta
  vec4 ktil = pait + pbt;
  vec4 k = pa - pi + pb;
  vec4 q1 = pz;
  vec4 q2 = (ktil + k) * (-2.) * (1. / (ktil + k).m2()) * ((ktil + k) * pz);
  vec4 q3 = k * (2.) * (1. / k.m2()) * (ktil * pz);
  pz = q1 + q2 + q3;

  // Define Mandelstam Variables (s_bar t_bar u_bar)
  // q(pa) qbar(pb) -> Z(pz) g(pi)
  // g(pa) qbar(pb) -> Z(pz) qbar(pi)
  // g(pa) q(pb) -> Z(pz) q(pi)
  double s_bar = (pa + pb).m2();
  double t_bar = (pb - pi).m2();
  double u_bar = (pa - pi).m2();

  // -------------------------------------------------------------------------
  // Calculate the cross-section

  // Calculate PDFs
  fl_a[idx] = !is_q2qg ? 21 : (ij_is_quark ? fl : -fl);
  fl_b[idx] = ij_is_quark ? -fl : fl;
  xa[idx] = ij_is_quark ? eta_q / x : eta_qbar / x;
  xb[idx] = ij_is_quark ? eta_qbar : eta_q;
  q2pdf[idx] = mu2;

  // Preamble
  double dxs_nlo;
  dxs_nlo = 2.;  // Two Possible Orientations (P(q) P(qbar) or P(qbar) P(q))
  dxs_nlo *= log((matrix->root_s * matrix->root_s) / pz.m2());  // from LO
  dxs_nlo *= (1. - x_min) * (1. - x);                           // x, v Sampling
  dxs_nlo *= 1. / (16. * M_PI);
  dxs_nlo *= 1. / pz.m2();  // Leftover
  double pd[5] = {0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279};
  dxs_nlo /= pd[ev.get_particle(0).get_pid() - 1];  // Flavour Selection
  dxs_nlo *= GeV_minus_2_to_pb;                     // units

  // Subtraction
  if (is_q2qg) {
    dxs_nlo *= 8. * M_PI * (*as)(mu2)*k_cf * me2_lo / pz.m2();
    dxs_nlo *= -2.;
  } else {
    dxs_nlo *= 8. * M_PI * (*as)(mu2)*k_tr * me2_lo / pz.m2();
    dxs_nlo *= (2. * pz.m2() - u_bar) / (s_bar);
    dxs_nlo *= 2.;  // two orientations
  }

  // 2 x types of event
  dxs_nlo *= 2.;

  // Adjust dxs for the nlo weighting
  if (matrix->ws != 0.) dxs_nlo /= matrix->ws;

  // -------------------------------------------------------------------------
  // Finalise the Event

  // Create the particle - can be gluon or quark!
  particle em;

  // set the momenta to the relevant partons
  // Need to use quark_is_element_0, ij_is_quark, and is_q2qg

  // quark emitter
  if (ij_is_quark) {
    // q->qg
    if (is_q2qg) {
      // quark
      ev.set_particle_mom(0, pa);
      ev.set_particle_acol(0, 2);
      ev.set_particle_eta(0, eta_q / x);
      // gluon
      em.set_pid(21);
      em.set_mom(pi);
      em.set_col(2);
      em.set_acol(1);
    }
    // q -> g qbar
    else {
      // gluon
      ev.set_particle_pid(0, 21);
      ev.set_particle_mom(0, pa);
      ev.set_particle_col(0, 2);
      ev.set_particle_acol(0, 1);
      ev.set_particle_eta(0, eta_q / x);
      // antiquark
      em.set_pid(-fl);
      em.set_mom(pi);
      em.set_acol(2);
    }
    // antiquark
    ev.set_particle_mom(1, pb);
  }
  // qbar emitter
  else {
    // qbar -> qbar g
    if (is_q2qg) {
      // antiquark
      ev.set_particle_mom(1, pa);
      ev.set_particle_col(1, 2);
      ev.set_particle_eta(1, eta_qbar / x);
      // gluon
      em.set_pid(21);
      em.set_mom(pi);
      em.set_col(1);
      em.set_acol(2);
    }
    // qbar -> g q
    else {
      // gluon
      ev.set_particle_pid(1, 21);
      ev.set_particle_mom(1, pa);
      ev.set_particle_col(1, 1);
      ev.set_particle_acol(1, 2);
      ev.set_particle_eta(1, eta_qbar / x);
      // quark
      em.set_pid(fl);
      em.set_mom(pi);
      em.set_col(2);
    }
    // quark
    ev.set_particle_mom(0, pb);
  }

  // Set the momenta of the Z Boson
  ev.set_particle_mom(2, pz);

  // Add the gluon to the event
  ev.set_particle(3, em);

  // increment emissions (important) !!!!!!
  ev.increment_emissions();

  if (isnan(dxs_nlo)) {
    dxs_nlo = 0.;
  }

  // Set the new cross-section
  ev.set_dxs(dxs_nlo);

  // -------------------------------------------------------------------------
  // Matching - Set the starting scale for the shower

  // This is t of the first emission!
  double t_in = (1 - x - v) * v / x * m.m2();

  ev.set_shower_t(t_in);
  ev.set_shower_c(2);

  return;
}

// -----------------------------------------------------------------------------
// S Event Generation - pp -> Z (Born + Virtual + Inseration + Collinear)

__global__ void c_terms(event *events, int n, matrix *matrix, int *fl_a,
                        int *fl_b, double *xa, double *xb, double *q2,
                        double *c_term) {
  /**
   * @brief Calculate the c term for each event
   *
   * @param events array of event records
   * @param n number of events
   * @param c_term array to store the c term for each event
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event &ev = events[idx];
  // ---------------------------------------------

  // If Real Correction, return (already done!)
  if (ev.get_size() == 4) {
    // Place with nonsense values that won't be used (Avoids PDF Error)
    fl_a[idx] = 1;
    fl_b[idx] = 1;
    xa[idx] = 0.5;
    xb[idx] = 0.5;
    q2[idx] = 1000.;
    return;
  }

  // ---------------------------------------------------------------------------
  // Need to do this again because new kernel

  // Get the matrix element squared
  double me2_lo = ev.get_me2();
  double dxs_lo = ev.get_dxs();

  // Momentum Fractions
  double eta_q = ev.get_particle(0).get_eta();
  double eta_qbar = ev.get_particle(1).get_eta();

  // Scale Choice = mz2, the only final state particle in Born!
  vec4 pz = ev.get_particle(2).get_mom();
  double mu2 = pz.m2();

  // flavour of the qqbar pair
  int fl = abs(ev.get_particle(0).get_pid());

  // ---------------------------------------------------------------------------
  // Now onto the C Term Generation

  // C Terms to be Calculated
  double c1 = 0.;   // End Points
  double c2g = 0.;  // PDF Ratios g
  double c2h = 0.;  // PDF Ratios h
  double c3 = 0.;   // Deltas

  // If our scale is not == pz.m2() this term is non-zero
  double log_mu = log(mu2 / pz.m2());

  // -------------------------------------------------------------------------
  // subtract end-points from the plus functions
  // i.e. - int{0 -> eta_i} (K + P) d eta', only the non-singular terms

  // q -> qg: EP = -2 (-pi^2/6 + (log(1 - eta_i))^2 + Li2(1 - eta_i))
  double ep_q2qg;
  ep_q2qg = log(1 - eta_q) * log_mu;
  ep_q2qg += (M_PI * M_PI) / 6;
  ep_q2qg -= pow(log(1. - eta_q), 2.);
  ep_q2qg -= dilogarithm(1. - eta_q);
  ep_q2qg *= 2.;

  double ep_qbar2qbarg;
  ep_qbar2qbarg = log(1 - eta_qbar) * log_mu;
  ep_qbar2qbarg += (M_PI * M_PI) / 6;
  ep_qbar2qbarg -= pow(log(1. - eta_qbar), 2.);
  ep_qbar2qbarg -= dilogarithm(1. - eta_qbar);
  ep_qbar2qbarg *= 2.;

  // No + terms for qbar -> g q or q -> g qbar
  // double ep_q2gqbar = 0.;
  // double ep_qbar2gq = 0.;

  // These are negative
  c1 = -(ep_q2qg + ep_qbar2qbarg);

  // -------------------------------------------------------------------------
  // Now, sample the intergration variable (eta') and pick one of the four
  // processes to evaluate the function for, and add them to the
  // cross-section

  // Decides which of the four integrals to do
  int r = static_cast<int>(ev.gen_random() * 4.) + 1;

  // Common vars
  double x, x_jac;
  if (r == 1 || r == 3) {
    x = eta_q + (1. - eta_q) * ev.gen_random();
    x_jac = 1. - eta_q;
  } else {
    x = eta_qbar + (1. - eta_qbar) * ev.gen_random();
    x_jac = 1. - eta_qbar;
  }

  // These both are identical really! Let's still leave it like this
  // to show there are two contributions
  double g, h;

  // q -> q g
  if (r == 1) {
    // Get the PDF ratio data
    fl_a[idx] = fl;
    fl_b[idx] = fl;
    xa[idx] = eta_q / x;
    xb[idx] = eta_q;
    q2[idx] = mu2;

    // f = g(x) * (pdf_ratio - 1) + h(x) * pdf_ratio
    g = 2. / (1. - x) * (log((1. - x) / x) + log(1. - x) - log_mu);
    h = -(1. + x) * (log((1. - x) / x) + log(1. - x) - log_mu);
    h += (1. - x);

    c2g = g * x_jac * 4.;
    c2h = h * x_jac * 4.;
  }

  // qbar -> qbar g
  else if (r == 2) {
    // Get the PDF ratio data
    fl_a[idx] = -fl;
    fl_b[idx] = -fl;
    xa[idx] = eta_qbar / x;
    xb[idx] = eta_qbar;
    q2[idx] = mu2;

    // f = g(x) * (pdf_ratio - 1) + h(x) * pdf_ratio
    g = 2. / (1. - x) * (log((1. - x) / x) + log(1. - x) - log_mu);
    h = -(1. + x) * (log((1. - x) / x) + log(1. - x) - log_mu);
    h += (1. - x);

    c2g = g * x_jac * 4.;
    c2h = h * x_jac * 4.;
  }

  // q -> g qbar
  else if (r == 3) {
    // Get the PDF ratio data
    fl_a[idx] = 21;
    fl_b[idx] = fl;
    xa[idx] = eta_q / x;
    xb[idx] = eta_q;
    q2[idx] = mu2;

    // f = g(x) * (pdf_ratio - 1) + h(x) * pdf_ratio
    g = 0.;
    h = 2. * x * (1. - x);
    h +=
        (1. - 2. * x + 2. * x * x) * (log((1. - x) / x) + log(1. - x) - log_mu);

    c2g = g * x_jac * 4. * k_tr / k_cf;
    c2h = h * x_jac * 4. * k_tr / k_cf;
  }

  // qbar -> g q
  else {
    // Get the PDF ratio data
    fl_a[idx] = 21;
    fl_b[idx] = -fl;
    xa[idx] = eta_qbar / x;
    xb[idx] = eta_qbar;
    q2[idx] = mu2;

    // f = g(x) * (pdf_ratio - 1) + h(x) * pdf_ratio
    g = 0.;
    h = 2. * x * (1. - x);
    h +=
        (1. - 2. * x + 2. * x * x) * (log((1. - x) / x) + log(1. - x) - log_mu);

    c2g = g * x_jac * 4. * k_tr / k_cf;
    c2h = h * x_jac * 4. * k_tr / k_cf;
  }

  // -------------------------------------------------------------------------
  // Finally, add the delta (1-x) terms

  // q -> q g: delta = -5. + 2pi^2/3
  double delta_q2qg = -5. + 2 * M_PI * M_PI / 3. - (3. / 2.) * log_mu;
  double delta_qbar2qbarg = -5. + 2 * M_PI * M_PI / 3. - (3. / 2.) * log_mu;

  // No + terms for qbar -> g q or q -> g qbar
  // double delta_q2gqbar = 0.;
  // double delta_qbar2gq = 0.;

  c3 = delta_q2qg + delta_qbar2qbarg;

  c_term[4 * idx + 0] = c1;
  c_term[4 * idx + 1] = c2g;
  c_term[4 * idx + 2] = c2h;
  c_term[4 * idx + 3] = c3;

  return;
}

__global__ void bvic_terms(event *events, int n, matrix *matrix, alpha_s *as,
                           double *c_term, double *xf_a, double *xf_b) {
  /**
   * @brief Combine the c term for each event
   *
   * @param events array of event records
   * @param n number of events
   * @param c_term array to store the c term for each event
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event &ev = events[idx];
  // ---------------------------------------------

  // If Real Correction, return (already done!)
  if (ev.get_size() == 4) {
    return;
  }

  // Calculate the Insertion Operator
  /**
   * V = - 2/e^2 - 3/e - 8 + pi^2
   * I(132) = 1/e^2 + 3/2e + 5 - pi^2/2
   * I(231) = 1/e^2 + 3/2e + 5 - pi^2/2
   * V + I(132) + I(231) = 2
   * Coeff = CF as(mu^2) / 2pi
   * Factor = Coeff * (V + I(132) + I(231)) = (CF as(mu^2) / 2pi) * 2
   *
   * dxs_nlo^(B + V + I) = dxs_lo * (1 + Factor)
   */
  double v_plus_i = 2.;

  // Get the c term
  double c1 = c_term[4 * idx + 0];
  double c2g = c_term[4 * idx + 1];
  double c2h = c_term[4 * idx + 2];
  double c3 = c_term[4 * idx + 3];

  // Calculate the total c term
  double pdf_ratio = xf_a[idx] / xf_b[idx];
  if (isnan(pdf_ratio) || isinf(pdf_ratio)) {
    pdf_ratio = 0.;
  }
  double c = c1 + c2g * (pdf_ratio - 1.) + c2h * pdf_ratio + c3;

  // Combine the B + V + I + C
  double mu2 = ev.get_particle(2).get_mom().m2();  // Scale Choice
  double coeff = k_cf * (*as)(mu2) / (2. * M_PI);
  double dxs_nlo = ev.get_dxs() * (1. + coeff * (v_plus_i + c));

  // Adjust dxs for the nlo weighting ws
  dxs_nlo /= (1. - matrix->ws);

  if (isnan(dxs_nlo)) {
    dxs_nlo = 0.;
  }

  if (isinf(dxs_nlo)) {
    dxs_nlo = 0.;
  }

  // Set the new cross-section
  ev.set_dxs(dxs_nlo);

  // -------------------------------------------------------------------------
  // Matching - Set the starting scale for the shower

  // Power Shower: Now allow emissions all the way up to root_s
  double sij = matrix->root_s * matrix->root_s;

  ev.set_shower_t(sij);
  ev.set_shower_c(1);
}

// -----------------------------------------------------------------------------
// Overall Wrapping Function

void lhc_nlo(thrust::device_vector<event> &dv_events, matrix *matrix,
             alpha_s *as, int blocks, int threads) {
  // Create PDF Evaluator
  pdf_wrapper pdf;

  // use a pointer to the device events
  event *d_events = thrust::raw_pointer_cast(dv_events.data());
  int n = dv_events.size();

  double *d_x_a;
  cudaMalloc(&d_x_a, n * sizeof(double));
  double *d_x_b;
  cudaMalloc(&d_x_b, n * sizeof(double));
  int *d_fl_a;
  cudaMalloc(&d_fl_a, n * sizeof(int));
  int *d_fl_b;
  cudaMalloc(&d_fl_b, n * sizeof(int));
  double *d_s_hat;
  cudaMalloc(&d_s_hat, n * sizeof(double));
  double *d_xf_a;
  cudaMalloc(&d_xf_a, n * sizeof(double));
  double *d_xf_b;
  cudaMalloc(&d_xf_b, n * sizeof(double));

  // Run LO Kernel
  debug_msg("running @lo_event");
  lo_event<<<blocks, threads>>>(d_events, n, matrix, d_fl_a, d_fl_b, d_x_a,
                                d_x_b, d_s_hat);
  sync_gpu_and_check("lo_event");

  // Launch the kernel
  debug_msg("running @h_event");
  h_event<<<blocks, threads>>>(d_events, n, matrix, as, d_fl_a, d_fl_b, d_x_a,
                               d_x_b, d_s_hat);
  sync_gpu_and_check("h_event");

  // Evalulate PDFs and multiply to dxs
  // - Real Emission would have changed the PDF params
  // - If no real Emission, PDF params stay the same as LO
  pdf.evaluate(d_fl_a, d_x_a, d_s_hat, d_xf_a, n, blocks, threads);
  pdf.evaluate(d_fl_b, d_x_b, d_s_hat, d_xf_b, n, blocks, threads);
  multiply_dxs_by_pdf<<<blocks, threads>>>(d_events, n, d_xf_a);
  multiply_dxs_by_pdf<<<blocks, threads>>>(d_events, n, d_xf_b);

  // Born + Virtual + Insertion + Collinear

  // Create 4xN array for c_terms
  double *d_c_term;
  cudaMalloc(&d_c_term, 4 * n * sizeof(double));
  cudaMemset(d_c_term, 0, 4 * n * sizeof(double));  // Initialize to zero

  debug_msg("running @c_terms");
  // Reuse the existing arrays instead of creating new _coll arrays
  c_terms<<<blocks, threads>>>(d_events, n, matrix, d_fl_a, d_fl_b, d_x_a,
                               d_x_b, d_s_hat, d_c_term);
  sync_gpu_and_check("c_terms");

  // Evaluate PDFs for collinear terms - reuse existing arrays
  pdf.evaluate(d_fl_a, d_x_a, d_s_hat, d_xf_a, n, blocks, threads);
  pdf.evaluate(d_fl_b, d_x_b, d_s_hat, d_xf_b, n, blocks, threads);

  // Do BVIC
  debug_msg("running @bvic_terms");
  bvic_terms<<<blocks, threads>>>(d_events, n, matrix, as, d_c_term, d_xf_a,
                                  d_xf_b);
  sync_gpu_and_check("bvic_terms");

  // Free only the c_term array (other arrays freed below)
  cudaFree(d_c_term);

  // Free the original arrays
  cudaFree(d_x_a);
  cudaFree(d_x_b);
  cudaFree(d_fl_a);
  cudaFree(d_fl_b);
  cudaFree(d_s_hat);
  cudaFree(d_xf_a);
  cudaFree(d_xf_b);
}