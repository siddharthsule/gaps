#include "decay_kinematics.cuh"

// ---------------------------------------------------------------------------
// Two-body throws

__device__ bool one_to_two_decay(vec4 p0, double m1, double m2, vec4& p1,
                                vec4& p2, double rho_1, double rho_2) {
  /**
   * @brief Two-body decay of p0 into m1 and m2, isotropic in its rest frame
   *
   * @param p0 momentum of the decaying particle
   * @param m1, m2 masses of the decay products
   * @param p1, p2 momenta of the decay products (output)
   * @param rho_1, rho_2 random numbers in [0, 1] for cos(theta) and phi
   * @return false where the decay is kinematically forbidden
   */

  // Mass of the decaying particle
  double m0 = p0.m();

  // Check if the decay is kinematically allowed (with some leeway)
  if (m0 < (m1 + m2) - 1e-9 * fmax(1., m1 + m2)) {
    printf("one_to_two_decay: m0 < m1 + m2!\n");
    return false;
  }

  // Momentum and energies of the decay products in the rest frame
  double p = p_star(m0, m1, m2);
  double e1 = (sqr(m0) + sqr(m1) - sqr(m2)) / (2. * m0);
  double e2 = (sqr(m0) + sqr(m2) - sqr(m1)) / (2. * m0);

  // Random isotropic directions for the decay products
  double ct = 2. * rho_1 - 1.;
  double st = sqrt(1. - ct * ct);
  double phi = 2. * M_PI * rho_2;

  // Unit vectors by construction: sin^2(theta) + cos^2(theta) = 1
  vec4 dir_p1(0, st * cos(phi), st * sin(phi), ct);
  vec4 dir_p2(0, -st * cos(phi), -st * sin(phi), -ct);

  p1 = vec4(e1, p * dir_p1[1], p * dir_p1[2], p * dir_p1[3]);
  p2 = vec4(e2, p * dir_p2[1], p * dir_p2[2], p * dir_p2[3]);

  // Boost the momenta back to the lab frame
  p1 = p0.boost_back(p1);
  p2 = p0.boost_back(p2);

  return true;
}

__device__ bool one_to_two_decay(vec4 p0, double m1, double m2, vec4& p1,
                                vec4& p2, vec4 axis) {
  /**
   * @brief Two-body decay of p0, the first product along axis in its rest frame
   *
   * @param p0 momentum of the decaying particle
   * @param m1, m2 masses of the decay products
   * @param p1, p2 momenta of the decay products (output)
   * @param axis rest-frame direction of p1, normalised here; +z if zero
   * @return false where the decay is kinematically forbidden
   */

  // Mass of the decaying particle
  double m0 = p0.m();

  // Kinematically allowed, to the same tolerance as the isotropic throw
  if (m0 < (m1 + m2) - 1e-9 * fmax(1., m1 + m2)) {
    printf("one_to_two_decay: m0 < m1 + m2!\n");
    return false;
  }

  // Momentum and energies of the decay products in the rest frame
  double p = p_star(m0, m1, m2);
  double e1 = (sqr(m0) + sqr(m1) - sqr(m2)) / (2. * m0);
  double e2 = (sqr(m0) + sqr(m2) - sqr(m1)) / (2. * m0);

  // Normalise the axis; fall back to +z if it is zero
  double axis_mag = axis.p();
  vec4 dir_p1 = (axis_mag > 0) ? vec4(0, axis[1] / axis_mag, axis[2] / axis_mag,
                                      axis[3] / axis_mag)
                               : vec4(0, 0, 0, 1);
  vec4 dir_p2 = vec4(0, -dir_p1[1], -dir_p1[2], -dir_p1[3]);

  p1 = vec4(e1, p * dir_p1[1], p * dir_p1[2], p * dir_p1[3]);
  p2 = vec4(e2, p * dir_p2[1], p * dir_p2[2], p * dir_p2[3]);

  // Boost the momenta back to the lab frame
  p1 = p0.boost_back(p1);
  p2 = p0.boost_back(p2);

  return true;
}

// ---------------------------------------------------------------------------
// n-body throw

__device__ bool one_to_n_decay(vec4 p0, const double* masses, int n,
                              vec4* momenta, event& ev) {
  /**
   * @brief n-body decay of p0, flat over the n-body phase space.
   *
   * Intermediate j is the system of the first n - 1 - j children.
   *
   * @param p0 momentum of the decaying particle
   * @param masses the n child masses
   * @param n the number of children, 2 to max_decay_products
   * @param momenta the n child momenta (output)
   * @param ev the event, for its random number stream
   * @return false where the decay is kinematically forbidden
   */

  // Mass of the decaying particle
  double m0 = p0.m();

  // Two children are one two-body throw, with no intermediate to draw
  if (n == 2) {
    if (m0 < masses[0] + masses[1]) {
      printf("one_to_n_decay: m0 < sum of child masses\n");
      return false;
    }

    double rho_1 = ev.gen_random();
    double rho_2 = ev.gen_random();
    return one_to_two_decay(p0, masses[0], masses[1], momenta[0], momenta[1],
                            rho_1, rho_2);
  }

  // Intermediate thresholds, summed inwards out so mins[j - 1] == mins[j] +
  // masses[n - 1 - j] exactly; an ulp of drift fails a decay at threshold
  double mins[max_decay_products];
  mins[n - 3] = masses[0] + masses[1];
  for (int j = n - 4; j >= 0; j--) {
    mins[j] = mins[j + 1] + masses[n - 2 - j];
  }

  // The sum of the child masses, to the same summation order
  double m_sum = mins[0] + masses[n - 1];

  // Check if the decay is kinematically allowed
  if (m0 < m_sum) {
    printf("one_to_n_decay: m0 < sum of child masses\n");
    return false;
  }

  // Mass the parent has left over, the range width every intermediate shares
  double a = m0 - m_sum;

  // The masses each intermediate reaches
  double maxs[max_decay_products];
  for (int j = 0; j < n - 2; j++) {
    maxs[j] = mins[j] + a;
  }

  // Bound on the weight, each momentum at its own endpoint
  double wt_max = p_star(maxs[n - 3], masses[0], masses[1]);
  for (int j = 0; j < n - 2; j++) {
    double outer = (j == 0) ? m0 : maxs[j - 1];
    wt_max *= p_star(outer, mins[j], masses[n - 1 - j]);
  }

  // Hit-or-miss for the intermediates, on the thresholds if the range closes
  double inter[max_decay_products];
  for (int j = 0; j < n - 2; j++) inter[j] = mins[j];

  if (wt_max > 0.) {
    bool accepted = false;

    for (int trial = 0; trial < max_phase_space_trials && !accepted; trial++) {
      // Offsets above the thresholds, drawn in a fixed order before use so the
      // CPU and GPU stay in sync
      double offsets[max_decay_products];
      for (int j = 0; j < n - 2; j++) offsets[j] = a * ev.gen_random();

      // Sort descending, so each intermediate fits inside the one outside it
      for (int j = 1; j < n - 2; j++) {
        double key = offsets[j];
        int k = j - 1;
        while (k >= 0 && offsets[k] < key) {
          offsets[k + 1] = offsets[k];
          k--;
        }
        offsets[k + 1] = key;
      }

      for (int j = 0; j < n - 2; j++) inter[j] = mins[j] + offsets[j];

      // The bound with each mass at the value drawn
      double wt = p_star(inter[n - 3], masses[0], masses[1]);
      for (int j = 0; j < n - 2; j++) {
        double outer = (j == 0) ? m0 : inter[j - 1];
        wt *= p_star(outer, inter[j], masses[n - 1 - j]);
      }

      // Accept-reject
      double rho = ev.gen_random();
      if (wt / wt_max > rho) accepted = true;
    }

    // Out of trials, the intermediates sit on their thresholds
    if (!accepted) {
      for (int j = 0; j < n - 2; j++) inter[j] = mins[j];
    }
  }

  // Walk the chain outermost first, each step splitting off one child
  vec4 p = p0;
  for (int j = 0; j < n - 2; j++) {
    vec4 next;
    double rho_1 = ev.gen_random();
    double rho_2 = ev.gen_random();
    if (!one_to_two_decay(p, inter[j], masses[n - 1 - j], next,
                          momenta[n - 1 - j], rho_1, rho_2)) {
      return false;
    }
    p = next;
  }

  // The innermost intermediate splits into the first two children
  double rho_1 = ev.gen_random();
  double rho_2 = ev.gen_random();
  return one_to_two_decay(p, masses[0], masses[1], momenta[0], momenta[1],
                          rho_1, rho_2);
}
