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
// LO Event Generation - pp -> Z
// NOTE: This kernel is no longer used directly. Instead, lhc_nlo now calls
// lhc_lo() wrapper which handles LO event generation and PDF evaluation.

// -----------------------------------------------------------------------------
// H Event Generation - pp -> Zj

__global__ void h_event(event* events, int n, matrix* matrix, alpha_s* as,
                        int* fl_a, int* fl_b, double* xa, double* xb,
                        double* q2pdf) {
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event& ev = events[idx];
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

  // Power Shower: Now allow emissions all the way up to root_s
  double sij = matrix->root_s * matrix->root_s;

  ev.set_shower_t(sij);
  ev.set_shower_c(2);

  return;
}

// -----------------------------------------------------------------------------
// S Event Generation - pp -> Z (Born + Virtual + Inseration + Collinear)

__global__ void c_terms(event* events, int n, matrix* matrix, int* fl_a,
                        int* fl_b, double* xa, double* xb, double* q2,
                        double* c_term) {
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
  event& ev = events[idx];
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
  // Now onto the Collinear Term Calculation

  /**
   * As in the paper, we break down the contributions into three components:
   * - A: delta(1-x) terms
   * - B(eta): plus function terms
   * - the g and h functions
   */

  // If our scale is not == pz.m2() this term is non-zero
  double log_mu = log(mu2 / pz.m2());

  // -------------------------------------------------------------------------
  // A Term Calculation

  // q -> q g and qbar -> qbar g
  double delta_q2qg = -5. + 2 * M_PI * M_PI / 3. - (3. / 2.) * log_mu;
  double delta_qbar2qbarg = -5. + 2 * M_PI * M_PI / 3. - (3. / 2.) * log_mu;

  // Define our A term
  double A = delta_q2qg + delta_qbar2qbarg;

  // -------------------------------------------------------------------------
  // B term Calculation

  // q -> qg and qbar -> qbar g
  double ep_q2qg = 0.;
  ep_q2qg += log(1 - eta_q) * log_mu;
  ep_q2qg += (M_PI * M_PI) / 6;
  ep_q2qg -= pow(log(1. - eta_q), 2.);
  ep_q2qg -= dilogarithm(1. - eta_q);
  ep_q2qg *= 2.;

  double ep_qbar2qbarg = 0.;
  ep_qbar2qbarg += log(1 - eta_qbar) * log_mu;
  ep_qbar2qbarg += (M_PI * M_PI) / 6;
  ep_qbar2qbarg -= pow(log(1. - eta_qbar), 2.);
  ep_qbar2qbarg -= dilogarithm(1. - eta_qbar);
  ep_qbar2qbarg *= 2.;

  // No + terms for qbar -> g q or q -> g qbar
  // double ep_q2gqbar = 0.;
  // double ep_qbar2gq = 0.;

  // Define our B term
  double B = ep_q2qg + ep_qbar2qbarg;

  // -------------------------------------------------------------------------
  // g and h functions Calculation
  // Now, sample the intergration variable (eta') and pick one of the four
  // processes to evaluate the function for, and add them to the
  // cross-section

  // Initialize the g and h contributions
  double g_contrib = 0.;
  double h_contrib = 0.;

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

    g_contrib = g * x_jac * 4.;
    h_contrib = h * x_jac * 4.;
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

    g_contrib = g * x_jac * 4.;
    h_contrib = h * x_jac * 4.;
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
    h += (x * x + (1. - x) * (1. - x)) *
         (log((1. - x) / x) + log(1. - x) - log_mu);

    g_contrib = g * x_jac * 4. * k_tr / k_cf;
    h_contrib = h * x_jac * 4. * k_tr / k_cf;
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
    h += (x * x + (1. - x) * (1. - x)) *
         (log((1. - x) / x) + log(1. - x) - log_mu);

    g_contrib = g * x_jac * 4. * k_tr / k_cf;
    h_contrib = h * x_jac * 4. * k_tr / k_cf;
  }

  // -------------------------------------------------------------------------
  // Store the c term components

  c_term[4 * idx + 0] = A;
  c_term[4 * idx + 1] = B;
  c_term[4 * idx + 2] = g_contrib;
  c_term[4 * idx + 3] = h_contrib;

  return;
}

__global__ void bvic_terms(event* events, int n, matrix* matrix, alpha_s* as,
                           double* c_term, double* xf_a, double* xf_b) {
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
  event& ev = events[idx];
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
  double A = c_term[4 * idx + 0];
  double B = c_term[4 * idx + 1];
  double g_contrib = c_term[4 * idx + 2];
  double h_contrib = c_term[4 * idx + 3];

  // Calculate the total c term
  double pdf_ratio = xf_a[idx] / xf_b[idx];
  if (isnan(pdf_ratio) || isinf(pdf_ratio)) {
    pdf_ratio = 0.;
  }
  double c = A - B + g_contrib * (pdf_ratio - 1.) + h_contrib * pdf_ratio;

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

void lhc_nlo(thrust::device_vector<event>& dv_events, matrix* matrix,
             alpha_s* as, pdf_wrapper* pdf, int blocks, int threads) {
              
  // use a pointer to the device events
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n = dv_events.size();

  double* d_x_a;
  cudaMalloc(&d_x_a, n * sizeof(double));
  double* d_x_b;
  cudaMalloc(&d_x_b, n * sizeof(double));
  int* d_fl_a;
  cudaMalloc(&d_fl_a, n * sizeof(int));
  int* d_fl_b;
  cudaMalloc(&d_fl_b, n * sizeof(int));
  double* d_s_hat;
  cudaMalloc(&d_s_hat, n * sizeof(double));
  double* d_xf_a;
  cudaMalloc(&d_xf_a, n * sizeof(double));
  double* d_xf_b;
  cudaMalloc(&d_xf_b, n * sizeof(double));

  // H Event Generation - Real Emission Corrections
  debug_msg("running @h_event");
  h_event<<<blocks, threads>>>(d_events, n, matrix, as, d_fl_a, d_fl_b, d_x_a,
                               d_x_b, d_s_hat);
  sync_gpu_and_check("h_event");

  // Evaluate PDFs for H events
  // - Real Emission would have changed the PDF params
  // - If no real Emission, PDF params stay the same as LO
  pdf->evaluate(d_fl_a, d_x_a, d_s_hat, d_xf_a, n, blocks, threads);
  pdf->evaluate(d_fl_b, d_x_b, d_s_hat, d_xf_b, n, blocks, threads);
  pdf->multiply_dxs_by_pdf(d_events, n, d_xf_a, blocks, threads);
  pdf->multiply_dxs_by_pdf(d_events, n, d_xf_b, blocks, threads);

  // Born + Virtual + Insertion + Collinear

  // Create 4xN array for c_terms
  double* d_c_term;
  cudaMalloc(&d_c_term, 4 * n * sizeof(double));
  cudaMemset(d_c_term, 0, 4 * n * sizeof(double));  // Initialize to zero

  debug_msg("running @c_terms");
  // Reuse the existing arrays instead of creating new _coll arrays
  c_terms<<<blocks, threads>>>(d_events, n, matrix, d_fl_a, d_fl_b, d_x_a,
                               d_x_b, d_s_hat, d_c_term);
  sync_gpu_and_check("c_terms");

  // Evaluate PDFs for collinear terms - reuse existing arrays
  pdf->evaluate(d_fl_a, d_x_a, d_s_hat, d_xf_a, n, blocks, threads);
  pdf->evaluate(d_fl_b, d_x_b, d_s_hat, d_xf_b, n, blocks, threads);

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