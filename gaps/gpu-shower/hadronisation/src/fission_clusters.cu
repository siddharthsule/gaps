#include "hadronisation.cuh"

__global__ void fission_pass(event* events, cluster_list* cls,
                             hadronisation* had, int n, int* n_active) {
  /**
   * @brief Perform ONE fission pass over every event.
   *
   * Threshold condition: M^clpow >= clmax^clpow + (m1 + m2)^clpow
   *
   * Follows the classic Herwig approach (arXiv:0803.0883): the cluster is
   * split collinearly along the quark direction in the cluster rest frame.
   *
   * The original code ran a `while` loop *inside* the kernel, so every thread
   * had to wait for the single hardest event in the batch to finish splitting.
   * Instead we do just ONE pass here and let the host repeat the launch, the
   * same bulk-synchronous style the shower uses (see run_hadronisation). Each
   * pass, one thread handles one event; an event that still has a cluster above
   * threshold adds itself to `n_active`, telling the host another pass is
   * needed. Once no event does so, hadronisation fission is complete.
   *
   * @param events   the events
   * @param cls      per-event cluster lists
   * @param had      hadronisation parameters
   * @param n        number of events
   * @param n_active atomic counter: number of events still splitting this pass
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

  // Start with all_below = true, set to false if we encounter a cluster that
  // needs fissioning. If all clusters are below threshold, this event is done.
  bool all_below = true;

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

    // Above threshold, so another pass will be needed
    all_below = false;

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
  // If any cluster was still above threshold this pass, flag the event as
  // active so the host runs another pass.

  if (!all_below) {
    atomicAdd(n_active, 1);
  }
}
