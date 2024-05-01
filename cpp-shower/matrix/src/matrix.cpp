#include "matrix.h"

// Constructor
Matrix::Matrix(double alphas, double ecms)
    : alphas(alphas),
      ecms(ecms),
      MZ2(pow(91.1876, 2.)),
      GZ2(pow(2.4952, 2.)),
      alpha(1. / 128.802),
      sin2tw(0.22293),
      amin(1.e-10),
      ye(0.5),
      ze(0.01),
      ws(0.25) {}

// Leading Order Matrix Element Generation
double Matrix::ME2(int fl, double s, double t) {
  /**
   * Matrix Element Squared
   * ----------------------
   *
   * Matrix element squared for massess 2x2 Scattering via a virtual photon
   * or a Z boson. See Peskin and Schroeder, Page 156.
   */
  double qe = -1.;
  double ae = -0.5;
  double ve = ae - 2. * qe * sin2tw;
  double qf = (fl == 2 || fl == 4) ? 2. / 3. : -1. / 3.;
  double af = (fl == 2 || fl == 4) ? 0.5 : -0.5;
  double vf = af - 2. * qf * sin2tw;
  double kappa = 1. / (4. * sin2tw * (1. - sin2tw));
  double chi1 = kappa * s * (s - MZ2) / (pow(s - MZ2, 2.) + GZ2 * MZ2);
  double chi2 = pow(kappa * s, 2.) / (pow(s - MZ2, 2.) + GZ2 * MZ2);
  double term1 = (1. + pow(1. + 2. * t / s, 2.)) *
                 (pow(qf * qe, 2.) + 2. * (qf * qe * vf * ve) * chi1 +
                  (ae * ae + ve * ve) * (af * af + vf * vf) * chi2);
  double term2 = (1. + 2. * t / s) * (4. * qe * qf * ae * af * chi1 +
                                      8. * ae * ve * af * vf * chi2);
  return pow(4. * M_PI * alpha, 2.) * 3. * (term1 + term2);
}

// Generate a point
void Matrix::GenerateLOPoint(Event &ev) {
  thread_local std::random_device rd;
  thread_local std::mt19937 gen(rd());

  // Same seed option. Turn off by commenting when not in use!
  // Having an if statement if no seed is given would not be a fair comparison
  // to the GPU, so commented out is better for now. Maybe in the future.
  // thread_local std::mt19937 gen(seed);

  std::uniform_real_distribution<> dis(0., 1.);  // Uniform distribution

  // Flavour
  int fl = std::rand() % 5 + 1;  // Faster than using dis(gen)

  Vec4 pa, pb, p1, p2;
  double p0;

  // ---------------------------------------------------------------------------
  // e+e- -> q qbar Kinematics

  // -------
  // Scalars
  p0 = 0.5 * ecms;
  double ct = 2. * dis(gen) - 1.;
  double st = std::sqrt(1. - ct * ct);
  double phi = 2. * M_PI * dis(gen);
  // -------
  // Vectors
  pa = Vec4(p0, 0., 0., p0);                                          // e+
  pb = Vec4(p0, 0., 0., -p0);                                         // e-
  p1 = Vec4(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);     // q
  p2 = Vec4(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);  // qbar

  // ---------------------------------------------------------------------------
  // DIS Kinematics - in Breit Frame

  // -------
  // Scalars
  p0 = 0.5 * ecms;
  double xi = dis(gen);        // Momentum Fraction of Quark
  double xip0 = xi * p0;       // xi = Momentum Fraction in DIS
  double q = dis(gen) * xip0;  // Momentum Transfer
  // -------
  // Vectors
  pa = Vec4(xip0, 0., 0., xip0);       // q (initial)
  pb = Vec4(xip0, 0., 0., -xip0);      // e- (initial)
  p1 = Vec4(xip0, 0., 0., -xip0 + q);  // q (final)
  p2 = Vec4(xip0, 0., 0., xip0 - q);   // e- (final)

  // ---------------------------------------------------------------------------
  // LHC Kinematics

  // -------
  // Scalars
  double eta1 = dis(gen);  // Momentum Fraction of Quark 1
  double eta2 = eta1;      // Momentum Fraction of Quark 2, temporarily the same
  double pB = 0.5 * ecms;  // Beam Centre of Mass
  p0 = 0.5 * ((eta1 * p0) + (eta2 * p0));  // Quark Centre of Mass
  // -------
  // Vectors
  pa = Vec4(eta1 * pB, 0., 0., eta1 * pB);   // q
  pb = Vec4(eta2 * pB, 0., 0., -eta2 * pB);  // qbar
  p1 = Vec4(p0, 0., 0., p0);                 // e+
  p2 = Vec4(p0, 0., 0., -p0);                // e-

  /**
   * Partonic Matrix Element for massless 2 -> 2 scattering
   * ------------------------------------------------------
   *
   * QCD and Collider Physics, Page 89: The Matrix Element Squared for
   * systems of 2 -> 2 Scatterinf are given by:
   *
   * |M|^2 ~ t^2 + u^2 / s^2 ~ (s channel)
   * |M|^2 ~ s^2 + u^2 / t^2 ~ (t channel)
   *
   * nb: s = (pa + pb)^2, t = (pa - p1)^2, u = (pa - p2)^2
   *
   * All in all, the matrix element squared is the same, bar the
   * kinematics we plug into the system. So we can use the code above
   * with the right choices of s, t and u for EE, DIS and LHC physics!
   */
  double lome = ME2(fl, (pa + pb).M2(), (pa - p1).M2());
  // double lome = ME2(fl, (pa - p1).M2(), (pa + pb).M2()); // DIS, s <-> t

  /**
   * Conversion to Diff Cross Section
   * ---------------------------------
   *
   * The conversion to the differential cross section is given by:
   *
   * dsigma/dt = SumAvg(ME2) / 16 pi s^2
   *
   * We also multiply by 5 to account for the 5 flavours, and 3.89379656e8 to
   * convert from GeV^-2 to pb.
   *
   * For e+e- and q qbar, s^2 = (pa + pb)^2 = Q^2
   * For DIS, s^2 = (pa + pb)^2 = ???
   */
  double dxs = lome / (16. * M_PI * std::pow(ecms, 2.));
  dxs *= 5. * 3.89379656e8;

  // Minus Momentum for initial state such that the sum is zero
  Parton p[4] = {Parton(-11, -pa, 0, 0), Parton(11, -pb, 0, 0),
                 Parton(fl, p1, 1, 0), Parton(-fl, p2, 0, 1)};

  // DIS
  // Parton p[4] = {Parton(fl, -pa, 1, 0), Parton(11, -pb, 0, 0),
  //                Parton(fl, p1, 1, 0), Parton(-fl, p2, 0, 1)};

  // NOTE: ColorConnected will need changing for Initial-Final Dipoles!

  // Set the Partons
  for (int i = 0; i < 4; i++) {
    ev.SetParton(i, p[i]);
  }

  // Set the ME Params
  ev.SetDxs(dxs);
  ev.SetHard(4);
}