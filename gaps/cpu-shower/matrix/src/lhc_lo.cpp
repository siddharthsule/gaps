#include "matrix.h"

void matrix::lhc_lo(event& ev) {
  /**
   * @brief generate the leading order interaction p p -> e+ e-
   *
   * The Leading order differential cross section for p p -> e+ e- is
   * given by:
   *
   * dsigma = drho1 drho2 drho3 drho4 * etaf(eta_a, s_hat) etaf(eta_b, s_hat)
   *          * (1/s_hat) * (I+ - I-)((s_hat - mz2)2 + mz2 gz2)/(mz gz)
   *          * ln(s/s_hat) * 1/(2 s_hat) 1/(8pi) * 2 * 1/N_c
   *          * 1/4 Sum_spins |M|^2
   *
   * @param ev the event object
   */

  // LHC: Bias Flavour based on contribution (based on cross section)
  // Flavours:    [down,       up,         strange,    charm,      bottom    ]
  // Percentages: [0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279]
  // Cumulatives: [0.36677679, 0.80186718, 0.91861688, 0.97129721, 1.00000000]
  double pd[5] = {0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279};
  double pc[5] = {0.36677679, 0.80186718, 0.91861688, 0.97129721, 1.00000000};

  // For unbiased Flavour Production
  // double pd[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
  // double pc[5] = {0.2, 0.4, 0.6, 0.8, 1.0};

  // Generate Flavor
  int fl;
  double r = ev.gen_random();
  for (int i = 0; i < 5; ++i) {
    if (r < pc[i]) {
      fl = i + 1;  // Flavors are 1-based
      break;
    }
  }

  // Four random numbers for s_hat, y, cos(theta) and phi
  double rho_1 = ev.gen_random();
  double rho_2 = ev.gen_random();
  double rho_3 = ev.gen_random();
  double rho_4 = ev.gen_random();

  // Generate s_hat
  double I_a = atan(((mz_cut_a * mz_cut_a) - mz2) / (mz * gz));
  double I_b = atan(((mz_cut_b * mz_cut_b) - mz2) / (mz * gz));
  double I = I_a + (I_b - I_a) * rho_1;
  double s_hat = (mz * gz * tan(I)) + mz2;

  // generate y
  double lim = min(100., 0.5 * log((root_s * root_s) / s_hat));
  double y = -lim + 2. * lim * rho_2;

  // Generate cos(theta), phi
  double ct = -1. + 2. * rho_3;
  double phi = 2. * M_PI * rho_4;

  // Kinematics
  // Momentum Fractions
  double eta_a = sqrt(s_hat) / root_s * exp(y);
  double eta_b = sqrt(s_hat) / root_s * exp(-y);
  // Quark / AntiQuark Momentum
  vec4 pa = vec4(eta_a * 0.5 * root_s, 0., 0., eta_a * 0.5 * root_s);
  vec4 pb = vec4(eta_b * 0.5 * root_s, 0., 0., -eta_b * 0.5 * root_s);
  // Z Momentum
  vec4 q = pa + pb;
  // Z Rest Frame
  vec4 q_rest = q.boost(q);
  // Electron Positron in Rest Frame
  double p0 = 0.5 * q_rest.m();
  double st = sqrt(1. - ct * ct);  // sin(theta)
  vec4 pr1 = vec4(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);
  vec4 pr2 = vec4(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);
  // Boost back to Lab Frame
  vec4 p1 = q.boost_back(pr1);
  vec4 p2 = q.boost_back(pr2);

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

  // Generate the particles
  particle p[4] = {particle(), particle(), particle(), particle()};
  p[0] = particle(fl, pa, 0, 1, eta_a);
  p[1] = particle(-fl, pb, 1, 0, eta_b);
  p[2] = particle(-11, p1, 0, 0);
  p[3] = particle(11, p2, 0, 0);

  // Calculate the Matrix Element
  double lome = me2(fl, (pa + pb).m2(), (pa - p1).m2());
  lome *= 1 / k_nc;  // Three possible initial colour states
  lome *= 2.;        // Two Possible Orientations

  // Evaluate PDFs
  double pdf_a = pdf->xfxQ2(fl, eta_a, s_hat);   // x_a f(x_a, s_hat)
  double pdf_b = pdf->xfxQ2(-fl, eta_b, s_hat);  // x_b f(x_b, s_hat)

  // Calculate the Cross Section
  double dxs;
  dxs = pdf_a * pdf_b;
  dxs *= (1. / s_hat);
  dxs *= (I_b - I_a);
  dxs *= (pow((s_hat - mz2), 2) + mz2 * gz2) / (mz * gz);
  dxs *= log((root_s * root_s) / s_hat);
  dxs *= (1. / (2. * s_hat)) * (1. / (8. * M_PI)) * lome;
  dxs *= GeV_minus_2_to_pb;  // units
  dxs /= pd[abs(fl) - 1];    // Flavour Selection

  // Set the particles
  for (int i = 0; i < 4; i++) {
    ev.set_particle(i, p[i]);
  }
  ev.set_hard(4);

  // Store the Matrix Element and Cross Section
  ev.set_me2(lome);
  ev.set_dxs(dxs);
}