#include "hadronisation.cuh"

__device__ void hadronisation::kallen(vec4 p0, double m1, double m2, vec4& p1,
                                      vec4& p2, double rho_1,
                                      double rho_2) const {
  /**
   * @brief Kallen function for decay kinematics
   *
   * Creates the kinematics for a two-body decay of a particle with
   * momentum p0 into two particles with masses m1 and m2. If axis is
   * provided, it is used as the direction of the first decay product.
   * Otherwise, isotropic decay is assumed.
   *
   * @param m0 mass of the decaying particle
   * @param m1 mass of the first decay product
   * @param m2 mass of the second decay product
   * @param p0 momentum of the decaying particle
   * @param p1 momentum of the first decay product (output)
   * @param p2 momentum of the second decay product (output)
   */

  // Mass of the decaying particle
  double m0 = p0.m();

  // Check if the decay is kinematically allowed
  if (m0 < m1 + m2) {
    printf("kallen: m0 < m1 + m2!\n");
    printf("m0 = %f, m1 = %f, m2 = %f\n", m0, m1, m2);
    return;
  }

  // Generate the Kallen function for the decay
  double lam = (sqr(m0) - sqr(m1 + m2)) * (sqr(m0) - sqr(m1 - m2));

  // Calculate the momentum of the decay products
  double p = 0.5 * sqrt(lam) / m0;

  // Calculate the energy of the decay products
  double e1 = (sqr(m0) + sqr(m1) - sqr(m2)) / (2. * m0);
  double e2 = (sqr(m0) + sqr(m2) - sqr(m1)) / (2. * m0);

  // Create random isotropic directions for the decay products
  double ct = 2. * rho_1 - 1.;
  double st = sqrt(1. - ct * ct);
  double phi = 2. * M_PI * rho_2;

  // Directions of the decay products (unit vectors by construction: sin²θ +
  // cos²θ = 1)
  vec4 dir_p1(0, st * cos(phi), st * sin(phi), ct);
  vec4 dir_p2(0, -st * cos(phi), -st * sin(phi), -ct);

  // Calculate the momentum of the decay products in the rest frame
  p1 = vec4(e1, p * dir_p1[1], p * dir_p1[2], p * dir_p1[3]);
  p2 = vec4(e2, p * dir_p2[1], p * dir_p2[2], p * dir_p2[3]);

  // Boost the momenta back to the lab frame
  p1 = p0.boost_back(p1);
  p2 = p0.boost_back(p2);
}

__device__ void hadronisation::kallen(vec4 p0, double m1, double m2, vec4& p1,
                                      vec4& p2, vec4 axis) const {
  /**
   * @brief Kallen function for decay kinematics with a fixed axis direction.
   *
   * Like the isotropic overload but instead of random angles the first decay
   * product is emitted along the provided unit three-vector `axis` (E=0) in
   * the rest frame of p0. Used for collinear cluster fission.
   *
   * @param p0   four-momentum of the decaying particle (lab frame)
   * @param m1   mass of the first decay product
   * @param m2   mass of the second decay product
   * @param p1   output: four-momentum of the first decay product
   * @param p2   output: four-momentum of the second decay product
   * @param axis unit three-vector (E=0) giving the decay direction in the
   *             rest frame of p0
   */

  // Mass of the decaying particle
  double m0 = p0.m();

  // Check if the decay is kinematically allowed
  if (m0 < m1 + m2) {
    printf("kallen: m0 < m1 + m2!\n");
    printf("m0 = %f, m1 = %f, m2 = %f\n", m0, m1, m2);
    return;
  }

  // Calculate the Kallen function for the decay
  double lam = (sqr(m0) - sqr(m1 + m2)) * (sqr(m0) - sqr(m1 - m2));

  // Calculate the momentum of the decay products
  double p = 0.5 * sqrt(lam) / m0;

  // Calculate the energy of the decay products
  double e1 = (sqr(m0) + sqr(m1) - sqr(m2)) / (2. * m0);
  double e2 = (sqr(m0) + sqr(m2) - sqr(m1)) / (2. * m0);

  // Normalise the axis; fall back to +z if it is zero
  double axis_mag = axis.p();
  vec4 dir_p1 = (axis_mag > 0) ? vec4(0, axis[1] / axis_mag, axis[2] / axis_mag,
                                      axis[3] / axis_mag)
                               : vec4(0, 0, 0, 1);
  vec4 dir_p2 = vec4(0, -dir_p1[1], -dir_p1[2], -dir_p1[3]);

  // Calculate the momentum of the decay products in the rest frame
  p1 = vec4(e1, p * dir_p1[1], p * dir_p1[2], p * dir_p1[3]);
  p2 = vec4(e2, p * dir_p2[1], p * dir_p2[2], p * dir_p2[3]);

  // Boost the momenta back to the lab frame
  p1 = p0.boost_back(p1);
  p2 = p0.boost_back(p2);
}
