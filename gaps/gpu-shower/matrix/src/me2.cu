#include "matrix.cuh"

__device__ double matrix::me2(int fl, double s, double t) const {
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