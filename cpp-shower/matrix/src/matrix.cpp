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
void Matrix::GenerateLOPoint(Event &ev, int seed) {
  thread_local std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(0., 1.);  // Uniform distribution
  std::uniform_int_distribution<> dis_int(1, 5);

  double ct = 2. * dis(gen) - 1.;
  double st = std::sqrt(1. - ct * ct);
  double phi = 2. * M_PI * dis(gen);

  int fl = dis_int(gen);

  Vec4 p1(1, st * std::cos(phi), st * std::sin(phi), ct);
  p1 = p1 * ecms / 2.;
  Vec4 p2(p1[0], -p1[1], -p1[2], -p1[3]);
  Vec4 pa(ecms / 2., 0., 0., ecms / 2.);
  Vec4 pb(ecms / 2., 0., 0., -ecms / 2.);

  double lome = ME2(fl, (pa + pb).M2(), (pa - p1).M2());
  double dxs =
      5. * lome * 3.89379656e8 / (8. * M_PI) / (2. * std::pow(ecms, 2.));

  Parton p[4] = {Parton(-11, -pa, 0, 0, true, false),
                 Parton(11, -pb, 0, 0, true, false), Parton(fl, p1, 1, 0),
                 Parton(-fl, p2, 0, 1)};

  // Set the Partons
  for (int i = 0; i < 4; i++) {
    ev.SetParton(i, p[i]);
  }

  // Set the ME Params
  ev.SetDxs(dxs);
  ev.SetHard(4);
  ev.SetInitial(2);
  ev.SetNonParton(2);
}