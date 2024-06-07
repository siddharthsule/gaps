#include "matrix.h"

// constructor
matrix::matrix(double alphas, double ecms)
    : alphas(alphas),
      ecms(ecms),
      mz2(pow(91.1876, 2.)),
      gz2(pow(2.4952, 2.)),
      alpha(1. / 128.802),
      sin2tw(0.22293),
      amin(1.e-10),
      ye(0.5),
      ze(0.01),
      ws(0.25) {}

// leading order matrix element generation
double matrix::me2(int fl, double s, double t) {
  /**
   * matrix element squared
   * ----------------------
   *
   * matrix element squared for massess 2x2 scattering via a virtual photon
   * or a z boson. see peskin and schroeder, page 156.
   */
  double qe = -1.;
  double ae = -0.5;
  double ve = ae - 2. * qe * sin2tw;
  double qf = (fl == 2 || fl == 4) ? 2. / 3. : -1. / 3.;
  double af = (fl == 2 || fl == 4) ? 0.5 : -0.5;
  double vf = af - 2. * qf * sin2tw;
  double kappa = 1. / (4. * sin2tw * (1. - sin2tw));
  double chi1 = kappa * s * (s - mz2) / (pow(s - mz2, 2.) + gz2 * mz2);
  double chi2 = pow(kappa * s, 2.) / (pow(s - mz2, 2.) + gz2 * mz2);
  double term1 = (1. + pow(1. + 2. * t / s, 2.)) *
                 (pow(qf * qe, 2.) + 2. * (qf * qe * vf * ve) * chi1 +
                  (ae * ae + ve * ve) * (af * af + vf * vf) * chi2);
  double term2 = (1. + 2. * t / s) * (4. * qe * qf * ae * af * chi1 +
                                      8. * ae * ve * af * vf * chi2);
  return pow(4. * M_PI * alpha, 2.) * 3. * (term1 + term2);
}

// generate a point
void matrix::generate_lo_point(event &ev) {
  // PRNG - Get at the start of Function, Set at the end
  unsigned long seed = ev.get_seed();
  double rand = ev.get_rand();

  // Flavour
  int fl = static_cast<int>(rand) % 5 + 1;
  update_rng(seed, rand);

  vec4 pa, pb, p1, p2;
  double p0;

  // ---------------------------------------------------------------------------
  // e+e- -> q qbar kinematics

  // -------
  // scalars
  p0 = 0.5 * ecms;
  double ct = 2. * rand - 1.;
  update_rng(seed, rand);
  double st = std::sqrt(1. - ct * ct);
  double phi = 2. * M_PI * rand;
  update_rng(seed, rand);

  // -------
  // vectors
  pa = vec4(p0, 0., 0., p0);                                          // e+
  pb = vec4(p0, 0., 0., -p0);                                         // e-
  p1 = vec4(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);     // q
  p2 = vec4(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);  // qbar

  /*
  // ---------------------------------------------------------------------------
  // dis kinematics - in breit frame

  // -------
  // scalars
  p0 = 0.5 * ecms;
  double xi = dis(gen);        // momentum fraction of quark
  double xip0 = xi * p0;       // xi = momentum fraction in dis
  double q = dis(gen) * xip0;  // momentum transfer
  // -------
  // vectors
  pa = vec4(xip0, 0., 0., xip0);       // q (initial)
  pb = vec4(xip0, 0., 0., -xip0);      // e- (initial)
  p1 = vec4(xip0, 0., 0., -xip0 + q);  // q (final)
  p2 = vec4(xip0, 0., 0., xip0 - q);   // e- (final)

  // ---------------------------------------------------------------------------
  // lhc kinematics

  // -------
  // scalars
  double eta1 = dis(gen);  // momentum fraction of quark 1
  double eta2 = eta1;      // momentum fraction of quark 2, temporarily the same
  double p_b = 0.5 * ecms;                 // beam centre of mass
  p0 = 0.5 * ((eta1 * p0) + (eta2 * p0));  // quark centre of mass
  // -------
  // vectors
  pa = vec4(eta1 * p_b, 0., 0., eta1 * p_b);   // q
  pb = vec4(eta2 * p_b, 0., 0., -eta2 * p_b);  // qbar
  p1 = vec4(p0, 0., 0., p0);                   // e+
  p2 = vec4(p0, 0., 0., -p0);                  // e-
  */

  /**
   * partonic matrix element for massless 2 -> 2 scattering
   * ------------------------------------------------------
   *
   * qcd and collider physics, page 89: the matrix element squared for
   * systems of 2 -> 2 scatterinf are given by:
   *
   * |m|^2 ~ t^2 + u^2 / s^2 ~ (s channel)
   * |m|^2 ~ s^2 + u^2 / t^2 ~ (t channel)
   *
   * nb: s = (pa + pb)^2, t = (pa - p1)^2, u = (pa - p2)^2
   *
   * all in all, the matrix element squared is the same, bar the
   * kinematics we plug into the system. so we can use the code above
   * with the right choices of s, t and u for ee, dis and lhc physics!
   */
  double lome = me2(fl, (pa + pb).m2(), (pa - p1).m2());
  // double lome = me2(fl, (pa - p1).m2(), (pa + pb).m2()); // dis, s <-> t

  /**
   * conversion to diff cross section
   * ---------------------------------
   *
   * the conversion to the differential cross section is given by:
   *
   * dsigma/dt = sum_avg(me2) / 16 pi s^2
   *
   * we also multiply by 5 to account for the 5 flavours, and 3.89379656e8 to
   * convert from ge_v^-2 to pb.
   *
   * for e+e- and q qbar, s^2 = (pa + pb)^2 = q^2
   * for dis, s^2 = (pa + pb)^2 = ???
   */
  double dxs = lome / (16. * M_PI * std::pow(ecms, 2.));
  dxs *= 5. * 3.89379656e8;

  // minus momentum for initial state such that the sum is zero
  parton p[4] = {parton(-11, -pa, 0, 0), parton(11, -pb, 0, 0),
                 parton(fl, p1, 1, 0), parton(-fl, p2, 0, 1)};

  // dis
  // parton p[4] = {parton(fl, -pa, 1, 0), parton(11, -pb, 0, 0),
  //                parton(fl, p1, 1, 0), parton(-fl, p2, 0, 1)};

  // note: color_connected will need changing for initial-final dipoles!

  // set the partons
  for (int i = 0; i < 4; i++) {
    ev.set_parton(i, p[i]);
  }

  // set the me params
  ev.set_dxs(dxs);
  ev.set_hard(4);

  // set the seed and rand back to the event
  ev.set_seed(seed);
  ev.set_rand(rand);
  printf("Seed: %lu, Rand: %f\n", seed, rand);
}