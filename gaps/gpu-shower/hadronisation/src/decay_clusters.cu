#include "hadronisation.cuh"

__global__ void decay_clusters(event* events, cluster_list* cls,
                               hadronisation* had, int n) {
  /**
   * @brief Decay clusters into hadrons.
   *
   * A quark-antiquark (or diquark-antidiquark) pair is drawn from the vacuum
   * to split each cluster into two hadrons via isotropic two-body kinematics.
   *
   * Flavour selection: quark or diquark drawn from the pwt[] pool
   *   (see select_qq_flavour for the weight-based selection)
   *
   * Meson path (quark selected):
   *   cluster(q, qbar) -> meson(q, -fl) + meson(fl, qbar)
   *   If kinematically forbidden, redraw a fresh flavour and retry.
   *
   * Baryon path (diquark selected):
   *   cluster(q, qbar) -> baryon(q, fl) + antibaryon(qbar, -fl)
   *   Falls back to same-flavour hadron + photon on lookup failure.
   *
   * Single-hadron fallback:
   *   The cluster decays directly to its same-flavour hadron plus one photon
   *   via exact two-body kinematics.
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

  int n_cl = cl.get_size();

  for (int i = 0; i < n_cl; i++) {
    cluster c = cl.get_cluster(i);
    int i1 = c.get_i1();  // quark index in the event record
    int i2 = c.get_i2();  // antiquark index in the event record

    int p1_pid = ev.get_particle(i1).get_pid();  // quark pid (positive)
    int p2_pid = ev.get_particle(i2).get_pid();  // antiquark pid (negative)

    vec4 cl_mom = c.get_mom();
    double M = cl_mom.m();

    // -------------------------------------------------------------------------
    // If not allowed to split with d-dbar, single hadron fallback

    int dummy_pid;
    double mass_pid1, mass_pid2;
    bool ok1_test = which_hadron(p1_pid, -1, dummy_pid, mass_pid1);
    bool ok2_test = ok1_test && which_hadron(1, p2_pid, dummy_pid, mass_pid2);

    // If not possible, single hadron fallback: cluster -> hadron + photon
    if (!ok2_test || M < mass_pid1 + mass_pid2) {
      int had_pid = 0;
      double had_mass = 0.;
      if (!which_hadron(p1_pid, p2_pid, had_pid, had_mass)) continue;

      // Convert directly to hadron + photon with exact two-body kinematics
      vec4 p_had, p_gamma;
      double rho1 = ev.gen_random();
      double rho2 = ev.gen_random();
      had->kallen(cl_mom, had_mass, 0.0, p_had, p_gamma, rho1, rho2);

      ev.set_particle(i1, particle(had_pid, p_had, 0, 0));
      ev.set_particle(i2, particle(22, p_gamma, 0, 0));
      continue;
    }

    // -------------------------------------------------------------------------
    // Find a valid hadron pair (meson or baryon) within max_attempts tries

    /**
     * For clusters with few open channels (e.g. charm/bottom, where the baryon
     * lookups always fail), most attempts here are wasted draws. This could be
     * sped up by enumerating the fixed vacuum pool (d/u/s + the 9 diquarks)
     * once per cluster flavour, caching which pairs pass which_hadron, and
     * sampling from that cached set with the pwt weights. The kinematic (M <
     * m1+m2) check still has to stay per-cluster. Left as a future upgrade —
     * for light clusters the lookup nearly always succeeds, so the current
     * trial-and-error is cheap.
     */

    const int max_attempts = 1000;
    for (int attempts = 0; attempts < max_attempts; ++attempts) {
      // Select Flavour (draw the two randoms in a fixed order so the CPU and
      // GPU stay in sync — argument evaluation order is unspecified in C++)
      double rho = ev.gen_random();
      double rho_diq = ev.gen_random();
      int fl = had->select_qq_flavour(rho, true, rho_diq);

      // Is Diquark? Must be all normal or anti quarks
      bool diq = is_diquark(fl);
      int m1_pid, m2_pid;
      double m1_mass, m2_mass;
      bool ok1, ok2;

      if (diq) {
        // Baryon: cluster(q, qbar) -> baryon(q, fl) + antibaryon(qbar, -fl)
        ok1 = which_hadron(p1_pid, fl, m1_pid, m1_mass);
        ok2 = ok1 && which_hadron(p2_pid, -fl, m2_pid, m2_mass);
      } else {
        // Meson: cluster(q, qbar) -> meson(q, -fl) + meson(fl, qbar)
        ok1 = which_hadron(p1_pid, -fl, m1_pid, m1_mass);
        ok2 = ok1 && which_hadron(fl, p2_pid, m2_pid, m2_mass);
      }

      // Skip if hadron lookup failed or cluster too light to decay
      if (!ok2 || M < m1_mass + m2_mass) continue;

      // Valid split found — do splitting
      vec4 h1, h2;
      double rho1 = ev.gen_random();
      double rho2 = ev.gen_random();
      had->kallen(cl_mom, m1_mass, m2_mass, h1, h2, rho1, rho2);

      // Update the event record with the new hadrons
      ev.set_particle(i1, particle(m1_pid, h1, 0, 0));
      ev.set_particle(i2, particle(m2_pid, h2, 0, 0));

      // Exit the attempts loop since we successfully decayed the cluster
      break;
    }
  }
}
