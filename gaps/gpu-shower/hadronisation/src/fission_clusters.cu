#include "hadronisation.cuh"

__global__ void fission_clusters(event* events, cluster_list* cls,
                                 hadronisation* had, int n) {
  /**
   * @brief Fission clusters that are too heavy into two lighter clusters.
   *
   * Threshold condition: M^clpow >= clmax^clpow + (m1 + m2)^clpow
   *
   * Follows the classic Herwig approach (arXiv:0803.0883): the cluster is
   * split collinearly along the quark direction in the cluster rest frame.
   * The loop repeats until all splittable clusters fall below the threshold.
   *
   * @param ev the event to process
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  cluster_list& cl = cls[idx];
  // ---------------------------------------------

  const int max_attempts = 1000;

  bool all_below = false;
  while (!all_below) {
    // Snapshot cluster count so newly added clusters are checked next pass
    int n_cl = cl.get_size();

    for (int i = 0; i < n_cl; i++) {
      cluster c = cl.get_cluster(i);
      int i1 = c.get_i1();
      int i2 = c.get_i2();

      double M = c.get_mom().m();
      double m1 = ev.get_particle(i1).get_mom().m();
      double m2 = ev.get_particle(i2).get_mom().m();

      // Flavour-dependent parameters (light=0, charm=1, bottom=2)
      int fl1 = abs(ev.get_particle(i1).get_pid());
      int fl2 = abs(ev.get_particle(i2).get_pid());
      int tier = (fl1 == 5 || fl2 == 5) ? 2 : (fl1 == 4 || fl2 == 4) ? 1 : 0;

      // Threshold: M^clpow >= clmax^clpow + (m1+m2)^clpow
      if (pow(M, had->clpow[tier]) < pow(had->clmax[tier], had->clpow[tier]) +
                                         pow(m1 + m2, had->clpow[tier])) {
        continue;
      }

      // -----------------------------------------------------------------------
      // Find a kinematically valid split (up to max_attempts tries)

      double M1 = 0., M2 = 0.;
      int fl = 0;
      bool valid = false;

      // Using the range formula, check if the cluster can atleast split with
      // the addition of the lightest pair (d dbar)
      double mq_min = had->const_mass[0];  // mass of d quark
      if (M < m1 + m2 + 2.0 * mq_min) {
        continue;
      }

      for (int attempts = 0; attempts < max_attempts; ++attempts) {
        // Random quark flavor from pool
        fl = had->select_qq_flavour(ev.gen_random());
        double mq = had->const_mass[fl - 1];

        // Kinematic condition: M > m1 + m2 + 2*mq
        if (M < m1 + m2 + 2.0 * mq) continue;

        // Sample M1 and M2
        double range = M - m1 - m2 - 2.0 * mq;
        double pp = had->psplit[tier];
        M1 = (m1 + mq) + range * pow(ev.gen_random(), 1.0 / pp);
        M2 = (m2 + mq) + range * pow(ev.gen_random(), 1.0 / pp);

        // Check if M1 + M2 ≤ M
        if (M1 + M2 <= M) {
          valid = true;
          break;
        }
      }

      if (!valid) {
        continue;
      }

      double mq = had->const_mass[fl - 1];
      vec4 c_lab = c.get_mom();

      // -----------------------------------------------------------------------
      // Splitting axis: quark direction in the cluster rest frame

      vec4 q_rest = c_lab.boost(ev.get_particle(i1).get_mom());
      double qp = q_rest.p();
      vec4 axis = (qp > 0)
                      ? vec4(0, q_rest[1] / qp, q_rest[2] / qp, q_rest[3] / qp)
                      : vec4(0, 0, 0, 1);

      // -----------------------------------------------------------------------
      // Two-body kinematics: parent cluster → C1 + C2

      vec4 c1_lab, c2_lab;
      had->kallen(c_lab, M1, M2, c1_lab, c2_lab, axis);

      // C1 → original quark (q1) + new antiquark
      vec4 c_q1_lab_new, c_qbar_lab;
      had->kallen(c1_lab, m1, mq, c_q1_lab_new, c_qbar_lab, axis);

      // C2 → new quark + original antiquark (q2)
      vec4 c_q_lab, c_qbar2_lab_new;
      had->kallen(c2_lab, mq, m2, c_q_lab, c_qbar2_lab_new, axis);

      // -----------------------------------------------------------------------
      // Update particles and clusters

      // Colour index shared by the original cluster constituents
      int col_orig = ev.get_particle(i1).get_col();

      // Update original quark momentum (colour indices unchanged)
      ev.set_particle_mom(i1, c_q1_lab_new);

      // Add new antiquark: colour-connected to original quark
      int new_aq_idx = ev.get_size();
      ev.add_emission(particle(-fl, c_qbar_lab, 0, col_orig));

      // Add new quark: colour-connected to original antiquark
      int new_q_idx = ev.get_size();
      ev.add_emission(particle(fl, c_q_lab, col_orig, 0));

      // Update original antiquark momentum (colour indices unchanged)
      ev.set_particle_mom(i2, c_qbar2_lab_new);

      // Replace cluster i → (original quark + new antiquark)
      cl.set_cluster(i, cluster(i1, new_aq_idx,
                                ev.get_particle(i1).get_mom() +
                                    ev.get_particle(new_aq_idx).get_mom()));

      // Append cluster (new quark + original antiquark)
      cl.add_cluster(
          new_q_idx, i2,
          ev.get_particle(new_q_idx).get_mom() + ev.get_particle(i2).get_mom());
    }

    // -------------------------------------------------------------------------
    // Check if all non-exception clusters are now below threshold

    all_below = true;
    for (int j = 0; j < cl.get_size(); j++) {
      cluster c = cl.get_cluster(j);
      double M = c.get_mom().m();
      double m1 = ev.get_particle(c.get_i1()).get_mom().m();
      double m2 = ev.get_particle(c.get_i2()).get_mom().m();

      int fl1 = abs(ev.get_particle(c.get_i1()).get_pid());
      int fl2 = abs(ev.get_particle(c.get_i2()).get_pid());
      int tier = (fl1 == 5 || fl2 == 5) ? 2 : (fl1 == 4 || fl2 == 4) ? 1 : 0;

      if (pow(M, had->clpow[tier]) >= pow(had->clmax[tier], had->clpow[tier]) +
                                          pow(m1 + m2, had->clpow[tier])) {
        all_below = false;
        break;
      }
    }
  }
}
