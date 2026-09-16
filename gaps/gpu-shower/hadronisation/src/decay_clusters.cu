#include "hadronisation.cuh"

// ---------------------------------------------------------------------------
/**
 * About this file
 * ---------------
 * The cluster decay needs a lot of supporting machinery, which makes this one
 * of the largest files in the codebase. It is arranged as:
 *
 * - Hadron weight formulae (the physics)
 * - Candidate infrastructure, and the table built from it
 * - The cluster decay itself
 *
 * The key idea is to build a database of cluster decay products out of the
 * hadron database. Building it once and reading it back skips the work the
 * decay would otherwise redo for every cluster.
 */

// ---------------------------------------------------------------------------
// The light diquarks the vacuum pops, in the order the decay walks them:
// f1 f2 0 nJ with f1 >= f2 and nJ = 1 (scalar) or 3 (vector), less the
// identical flavours at nJ = 1 that Pauli forbids.

__device__ const int n_light_diquarks = 9;
__device__ const int light_diquarks[n_light_diquarks] = {1103, 2101, 2103, 2203,
                                              3101, 3103, 3201, 3203, 3303};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// HADRON WEIGHT
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Weight of a hadron from (q, qbar) or (q, diquark): 2J + 1 x meson mixing x
// baryon multiplet, zero if forbidden. The caller applies pwt and p_star.

__device__ double hadron_weight(int pid1, int pid2, int had_pid) {
  /**
   * @brief the weight a hadron carries in its pair's list
   *
   * 2J + 1, times the meson mixing or the baryon multiplet factor. The caller
   * applies pwt and p_star on top.
   *
   * @param pid1 the quark pid, either sign
   * @param pid2 the antiquark or diquark pid it pairs with
   * @param had_pid the hadron id, signed
   * @return the weight, zero where the pair cannot form the state
   */

  int a = abs(had_pid);
  int n1 = (a / 1000) % 10;
  int n2 = (a / 100) % 10;
  int n3 = (a / 10) % 10;
  int nj = a % 10;

  // Baryons
  if (n1 != 0) {
    // No excited states, no spin above 3/2, no Lambda-type ordering at
    // J = 3/2, no flavours out of descending order
    if (a > 9999 || nj > 4) return 0.;
    if (nj == 4 && n2 < n3) return 0.;
    if (n1 < n2 || n1 < n3) return 0.;

    // Spin weight 2J + 1, and Pauli: three identical flavours have no J = 1/2
    // partner, so the state that survives carries the weight of both
    return (n1 == n2 && n2 == n3) ? 1.5 * nj : nj;
  }

  // Mesons: the eta_c and the eta_b are not formed, being outside the mixing
  if (a == 441 || a == 551 || a == 100441 || a == 100551) return 0.;

  // Spin weight 2J + 1
  double weight = nj;

  // pi0 and higher: u ubar and d dbar share the state, so halve it
  if (n2 == n3 && n2 == 1) weight *= 0.5;

  // Mixing angles [degrees] measured from the light state, for the
  // pseudoscalar, vector, tensor and spin-3 nonets
  const int n_angles = 4;
  const double angles[n_angles] = {29.74, -1.74, 7.26, 4.46};

  // uds systems
  if (n2 == n3 && (n2 == 2 || n2 == 3)) {
    // A nonet carrying an excitation digit mixes ideally, as does one above
    // the tabulated angles, leaving the light state pure
    int index = (nj - 1) / 2;
    double phi =
        (a > 999 || index >= n_angles) ? 0. : angles[index] * M_PI / 180.;

    bool light_pair = (pid1 == 1 && pid2 == -1) || (pid1 == 2 && pid2 == -2);
    bool strange_pair = (pid1 == 3 && pid2 == -3);

    if (n2 == 2) {
      // Type 1, the 22nJ id (eta, omega, ...)
      if (light_pair)
        weight *= 0.5 * sqr(cos(phi));
      else if (strange_pair)
        weight *= sqr(sin(phi));
    } else {
      // Type 2, the 33nJ id (eta', phi, ...)
      if (light_pair)
        weight *= 0.5 * sqr(sin(phi));
      else if (strange_pair)
        weight *= sqr(cos(phi));
    }
  }

  return weight;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// CANDIDATE INFRASTRUCTURE
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

/**
 * Candidate structure
 * -------------------
 * Pre-computing every hadron the cluster decay can create drops the number of
 * weight computations sharply. An entry keeps the hadron's pid and mass, and
 * its weight: the spin nJ times the mixing and multiplet factors.
 */

struct candidate {
  int pid;        // signed
  double mass;    // GeV
  double weight;  // spin x mixing x multiplet
};

const int max_candidates = 24;

// ---------------------------------------------------------------------------
// Writing out one hadron of a list.

__device__ void add_candidate(int pid1, int pid2, int had_pid, candidate* out,
                          int& n) {
  /**
   * @brief Write out a hadron, where the table holds it and it carries weight.
   *
   * @param pid1 the quark pid, either sign
   * @param pid2 the antiquark or diquark pid it pairs with
   * @param had_pid the hadron id, signed; the table is keyed on its magnitude
   * @param out the candidates, max_candidates long (output)
   * @param n how many have been written (in and out)
   */

  double mass = hadron_mass(abs(had_pid));
  if (mass < 0.) return;

  double weight = hadron_weight(pid1, pid2, had_pid);
  if (weight <= 0.) return;

  if (n >= max_candidates) {
    printf("add_candidate: list full for %d %d\n", pid1, pid2);
    return;
  }

  out[n].pid = had_pid;
  out[n].mass = mass;
  out[n].weight = weight;
  n++;
}

// ---------------------------------------------------------------------------
// The hadrons a flavour pair forms; a pid2 in the diquark pool gives baryons.

__device__ int build_candidates(int pid1, int pid2, candidate* out) {
  /**
   * @brief Every hadron the two ids form with a mass and a non-zero weight.
   *
   * @param pid1 the quark pid, either sign
   * @param pid2 the antiquark or diquark pid it pairs with
   * @param out the candidates, max_candidates long (output)
   * @return the number of candidates written
   */

  int n = 0;

  // Mesons: a quark against an antiquark

  bool pid2_is_diquark = abs(pid2) >= 1000 && (abs(pid2) / 10) % 10 == 0;
  if (!pid2_is_diquark) {
    // The quark and the antiquark of the pair, either way round
    int q = (pid1 > 0) ? abs(pid1) : abs(pid2);
    int qbar = (pid1 > 0) ? abs(pid2) : abs(pid1);

    // The sign the two flavours give the id: a meson is a particle if its
    // heavier quark is an up-type quark or a down-type antiquark
    int sign;
    if (q == qbar) {
      sign = 1;
    } else if (q > qbar) {
      sign = (q % 2 == 0) ? 1 : -1;
    } else {
      sign = (qbar % 2 == 0) ? -1 : 1;
    }

    // The flavour digits, larger first. A diagonal pair runs over every nonet
    // its content reaches, so u ubar carries 11nJ, 22nJ and 33nJ
    int digits[3][2];
    int n_digits = 0;

    if (q == qbar) {
      if (q > 3) {
        // Charm and bottom keep their own digits
        digits[n_digits][0] = q;
        digits[n_digits][1] = q;
        n_digits++;
      } else {
        // Strange reaches the eta and eta' types, light also the pi type
        for (int d = (q == 3) ? 2 : 1; d <= 3; d++) {
          digits[n_digits][0] = d;
          digits[n_digits][1] = d;
          n_digits++;
        }
      }
    } else {
      // Two different flavours give just the one digit pair, larger first
      digits[n_digits][0] = max(q, qbar);
      digits[n_digits][1] = min(q, qbar);
      n_digits++;
    }

    // Every combination of nJ, nL and nr
    for (int d = 0; d < n_digits; d++) {
      int n2 = digits[d][0];
      int n3 = digits[d][1];

      for (int nj = 1; nj <= 9; nj += 2) {
        for (int nl = 0; nl <= 3; nl++) {
          for (int nr = 0; nr <= 3; nr++) {
            int pid = 100000 * nr + 10000 * nl + 100 * n2 + 10 * n3 + nj;
            add_candidate(pid1, pid2, sign * pid, out, n);
          }
        }
      }
    }

    return n;
  }

  // Baryons: a quark against one of the light diquarks

  int diq = abs(pid2);

  // A quark with a diquark makes the baryon, an antiquark with an
  // antidiquark its conjugate
  int sign = (pid1 > 0) ? 1 : -1;

  // The three flavours, largest first. The diquark already holds its two in
  // descending order, so only the quark needs placing among them
  int f1 = (diq / 1000) % 10;
  int f2 = (diq / 100) % 10;
  int q = abs(pid1);

  int n1, n2, n3;
  if (q >= f1) {
    n1 = q;
    n2 = f1;
    n3 = f2;
  } else if (q >= f2) {
    n1 = f1;
    n2 = q;
    n3 = f2;
  } else {
    n1 = f1;
    n2 = f2;
    n3 = q;
  }

  // The baryon id bar its nJ digit, which the couplings below supply
  int base = 1000 * n1 + 100 * n2 + 10 * n3;

  // Only two values of nJ are reachable: 2 for spin 1/2, 4 for spin 3/2
  if (diq % 10 == 1) {
    // A scalar diquark couples with the quark to J = 1/2 alone. Three distinct
    // flavours give the Lambda-type state, whose nq2 and nq3 the PDG swaps
    bool lambda_type = (n1 != n2 && n2 != n3);
    int pid = lambda_type ? 1000 * n1 + 100 * n3 + 10 * n2 + 2 : base + 2;
    add_candidate(pid1, pid2, sign * pid, out, n);
  } else if (n1 == n2 && n2 == n3) {
    // Pauli: three identical flavours have no J = 1/2 state
    add_candidate(pid1, pid2, sign * (base + 4), out, n);
  } else {
    // A vector diquark couples to both J = 1/2 and J = 3/2
    add_candidate(pid1, pid2, sign * (base + 2), out, n);
    add_candidate(pid1, pid2, sign * (base + 4), out, n);
  }

  return n;
}

// ---------------------------------------------------------------------------
// From the database, precompute the possible cluster decay products for baryons
// and mesons

// Our database of "candidates". A row is keyed on flavour counting from zero,
// so a quark pid q sits at q - 1 and an antiquark pid qbar at -qbar - 1; a
// diquark is named by its slot in light_diquarks.
struct candidate_table {
  // Mesons: 5 flavours x 5 flavours x the variants (nr, nL, nJ) reach
  candidate meson[5][5][max_candidates];

  // Baryons: 5 flavours x 9 diquarks x 2 variants (J = 1/2 or 3/2)
  candidate baryon[5][9][2];

  // Counters
  int n_meson[5][5];
  int n_baryon[5][9];
};

// hadron_mass reads a __device__ table, so the lists cannot be built on the
// host and pushed across: they are built on the device by the kernel below,
// which run_hadronisation launches once before the decay. The table is file
// local, so it stays out of the hadronisation object and out of its H2D copy.
__device__ candidate_table candidates;

__global__ void build_candidate_table() {
  /**
   * @brief fill the candidate table, one thread per flavour pair
   *
   * Launched once per run, before any kernel that reads the table.
   */

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = gridDim.x * blockDim.x;

  // The 25 meson pairs come first, then the 45 quark-diquark pairs
  int n_meson_pairs = 5 * 5;
  int n_pairs = n_meson_pairs + 5 * n_light_diquarks;

  for (int i = idx; i < n_pairs; i += stride) {
    // Mesons: the quark against each antiquark, straight into its row
    if (i < n_meson_pairs) {
      int q = i / 5;
      int qbar = i % 5;

      candidates.n_meson[q][qbar] =
          build_candidates(q + 1, -(qbar + 1), candidates.meson[q][qbar]);
      continue;
    }

    // Baryons: the quark against each popped diquark. A row holds only the
    // J = 1/2 and J = 3/2 states, so the list is built to one side first
    int j = i - n_meson_pairs;
    int q = j / n_light_diquarks;
    int d = j % n_light_diquarks;

    candidate rows[max_candidates];
    int n = build_candidates(q + 1, light_diquarks[d], rows);

    if (n > 2) {
      printf("build_candidate_table: %d %d forms %d baryons\n", q + 1,
             light_diquarks[d], n);
      n = 2;
    }

    candidates.n_baryon[q][d] = n;
    for (int k = 0; k < n; k++) candidates.baryon[q][d][k] = rows[k];
  }
}

// ---------------------------------------------------------------------------
// The lightest states a cluster reaches, read off the table

__device__ candidate lightest_meson(int q, int qbar) {
  /**
   * @brief the lightest meson a quark and an antiquark form
   *
   * @param q the quark pid, positive
   * @param qbar the antiquark pid, negative
   * @return the lightest candidate
   */

  const candidate* rows = candidates.meson[q - 1][-qbar - 1];
  int n = candidates.n_meson[q - 1][-qbar - 1];

  int best = 0;
  for (int i = 1; i < n; i++) {
    if (rows[i].mass < rows[best].mass) best = i;
  }

  return rows[best];
}

__device__ double lightest_meson_pair_mass(int q, int qbar) {
  /**
   * @brief the lightest meson pair mass a cluster reaches, over d, u, s pops
   *
   * @param q the cluster's quark pid, positive
   * @param qbar the cluster's antiquark pid, negative
   * @return the smallest sum
   */

  double best = INFINITY;

  for (int fl = 1; fl <= 3; fl++) {
    // cluster(q, qbar) -> meson(q, -fl) + meson(fl, qbar)
    double m = lightest_meson(q, -fl).mass + lightest_meson(fl, qbar).mass;

    if (m < best) best = m;
  }

  return best;
}

__device__ bool can_produce_baryon_pair(int q, int qbar, double m_cluster) {
  /**
   * @brief whether a cluster reaches a baryon pair at all
   *
   * Walks the pairs the popped diquarks open and stops at the first the
   * cluster is heavy enough to reach.
   *
   * @param q the cluster's quark pid, positive
   * @param qbar the cluster's antiquark pid, negative
   * @param m_cluster the cluster mass
   * @return true if any pop opens a pair the cluster reaches
   */

  for (int d = 0; d < n_light_diquarks; d++) {
    // cluster(q, qbar) -> baryon(q, diq) + antibaryon(qbar, -diq), the
    // antibaryon being the conjugate of the baryon the antiquark builds. A pop
    // the table holds no baryon for leaves a loop that does not run
    const candidate* rows1 = candidates.baryon[q - 1][d];
    int n1 = candidates.n_baryon[q - 1][d];
    const candidate* rows2 = candidates.baryon[-qbar - 1][d];
    int n2 = candidates.n_baryon[-qbar - 1][d];

    for (int a = 0; a < n1; a++) {
      for (int b = 0; b < n2; b++) {
        // Strictly below the cluster mass: a pair exactly at threshold has
        // p_star == 0, so it carries no weight and is never selected
        if (rows1[a].mass + rows2[b].mass < m_cluster) return true;
      }
    }
  }

  return false;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// MAIN HADRONISATION CODE
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Fallback Reshuffling

__device__ void fallback_reshuffling(event& ev, cluster_list& cl) {
  /**
   * @brief Reshuffle the momenta onto the masses the clusters decay with.
   *
   * A cluster below its lightest hadron pair takes its lightest hadron's mass;
   * k solves Sum( sqrt(k^2*|p|_i^2 + M_i^2) ) = Ecms.
   *
   * @param ev the event to reshuffle
   * @param cl the clusters formed from it
   */
  // Check for overflow
  if (ev.get_overflowed()) return;

  // Number of clusters
  int n_cl = cl.get_size();

  // Mark the particles that belong to a cluster
  bool in_cluster[max_particles] = {false};

  for (int i = 0; i < n_cl; i++) {
    in_cluster[cl.get_cluster(i).get_i1()] = true;
    in_cluster[cl.get_cluster(i).get_i2()] = true;
  }

  // The entries, and the mass each one is put on. The clusters come first,
  // in their own order, and the particles outside them follow
  vec4 momenta[max_particles];
  double masses[max_particles];
  int particle_index[max_particles];
  int n = 0;

  // Some events need no reshuffling at all
  bool needs_reshuffling = false;

  for (int i = 0; i < n_cl; i++) {
    cluster c = cl.get_cluster(i);
    double M = c.get_mom().m();

    int p1_pid = ev.get_particle(c.get_i1()).get_pid();
    int p2_pid = ev.get_particle(c.get_i2()).get_pid();

    momenta[n] = c.get_mom();

    // A cluster above the lightest pair it reaches keeps its own mass, and
    // one below it goes onto the lightest hadron its own flavours form
    if (M > lightest_meson_pair_mass(p1_pid, p2_pid)) {
      masses[n] = M;
    } else {
      masses[n] = lightest_meson(p1_pid, p2_pid).mass;
      needs_reshuffling = true;
    }

    n++;
  }

  // Early exit if no reshuffling needed
  if (!needs_reshuffling) return;

  // Anything outside a cluster keeps its own mass
  for (int i = 2; i < ev.get_size(); i++) {
    if (in_cluster[i]) continue;
    momenta[n] = ev.get_particle(i).get_mom();
    masses[n] = momenta[n].m();
    particle_index[n] = i;
    n++;
  }

  // Ecms is the total energy; the summed |p| sets the top of the bracket
  vec4 total;
  double p_sum = 0.;
  for (int i = 0; i < n; i++) {
    total = total + momenta[i];
    p_sum += momenta[i].p();
  }
  double ecms = total[0];

  // Every entry at rest leaves no momentum to scale
  if (p_sum == 0.) return;

  // No k exists once the masses alone exceed Ecms
  double mass_sum = 0.;
  for (int i = 0; i < n; i++) mass_sum += masses[i];

  if (mass_sum > ecms) {
    printf(
        "fallback_reshuffling: the masses the clusters decay with sum above "
        "Ecms, the momenta are left as they are\n");
    return;
  }

  // Bisection method to find k
  double klo = 0.0;
  double khi = fmax(1.0, ecms / p_sum);
  for (int iter = 0; iter < 50; iter++) {
    double kmid = 0.5 * (klo + khi);

    double f = 0.;
    for (int i = 0; i < n; i++) {
      f += sqrt(kmid * kmid * momenta[i].p2() + masses[i] * masses[i]);
    }
    f -= ecms;

    if (f > 0.0) {
      khi = kmid;
    } else {
      klo = kmid;
    }
    if (khi - klo < 1e-12) break;
  }
  double k = 0.5 * (klo + khi);

  // Reshuffle the momentum of every entry
  for (int i = 0; i < n; i++) {
    double e = sqrt(k * k * momenta[i].p2() + masses[i] * masses[i]);
    vec4 reshuffled(e, k * momenta[i][1], k * momenta[i][2], k * momenta[i][3]);

    if (i < n_cl) {
      cl.override_cluster_mom(i, reshuffled, ev);
    } else {
      ev.set_particle_mom(particle_index[i], reshuffled);
    }
  }
}

// ---------------------------------------------------------------------------
// The cluster decay itself: total the pairs one category opens, draw one of
// them, and decay the cluster into it.

__device__ double choose_decay(const double* pwt, int p1_pid, int p2_pid,
                               double m_cluster, bool want_baryon,
                               double target, int& h1_pid, double& h1_mass,
                               int& h2_pid, double& h2_mass) {
  /**
   * @brief Walk every hadron pair one category opens, summing their weights.
   *
   * W = w_pop * w_a * w_b * p_star(M, m_a, m_b), w_pop = pwt[q] for a quark
   * and pwt[3] * pwt[q] * pwt[q'] for a diquark. Every call walks the same
   * pairs in the same order, so the total and select walks agree bitwise.
   *
   * @param pwt the flavour weights, d u s and the diquark suppression
   * @param p1_pid the cluster's quark pid, positive
   * @param p2_pid the cluster's antiquark pid, negative
   * @param m_cluster the cluster mass
   * @param want_baryon pop a diquark rather than a quark
   * @param target the running sum at which to write out a pair; negative sums
   * @param h1_pid the first hadron's pid (output)
   * @param h1_mass the first hadron's mass (output)
   * @param h2_pid the second hadron's pid (output)
   * @param h2_mass the second hadron's mass (output)
   * @return the total weight walked
   */

  double total = 0.;

  // The cluster's own flavours, as the table is keyed on them
  int q = p1_pid - 1;
  int qbar = -p2_pid - 1;

  // The hadrons each side forms, read off the table for every flavour popped
  const candidate* list1;
  const candidate* list2;

  // What the vacuum pops: one of the light diquarks for a baryon pair, or a
  // d, u or s quark for a meson pair. light_diquarks is listed in the order
  // this walk visits them, with the f1 == f2 at nJ = 1 that Pauli forbids
  // already left out, so the pops need no guard of their own.
  int n_pop = want_baryon ? n_light_diquarks : 3;

  for (int p = 0; p < n_pop; p++) {
    double w_pop;
    int n1 = 0;
    int n2 = 0;
    int sign2 = 1;

    // Baryons
    if (want_baryon) {
      int diq = light_diquarks[p];
      int f1 = (diq / 1000) % 10;
      int f2 = (diq / 100) % 10;
      w_pop = pwt[3] * pwt[f1 - 1] * pwt[f2 - 1];

      // cluster(q, qbar) -> baryon(q, diq) + antibaryon(qbar, -diq), the
      // antibaryon being the conjugate of the baryon built here
      list1 = candidates.baryon[q][p];
      n1 = candidates.n_baryon[q][p];
      list2 = candidates.baryon[qbar][p];
      n2 = candidates.n_baryon[qbar][p];
      sign2 = -1;
    }

    // Mesons
    else {
      int f1 = p + 1;
      w_pop = pwt[f1 - 1];

      // cluster(q, qbar) -> meson(q, -f1) + meson(f1, qbar)
      list1 = candidates.meson[q][f1 - 1];
      n1 = candidates.n_meson[q][f1 - 1];
      list2 = candidates.meson[f1 - 1][qbar];
      n2 = candidates.n_meson[f1 - 1][qbar];
    }

    // For candidates pick using target = random(0,1)
    for (int a = 0; a < n1; a++) {
      for (int b = 0; b < n2; b++) {
        double m1 = list1[a].mass;
        double m2 = list2[b].mass;

        // Skip a pair the cluster is too light to reach
        if (m_cluster < m1 + m2) continue;

        // Generate Weight
        double w = w_pop * list1[a].weight * list2[b].weight *
                   p_star(m_cluster, m1, m2);

        // A pair of zero weight, one exactly at threshold, is never
        // selected
        if (w <= 0.) continue;

        total += w;
        if (target >= 0. && total >= target) {
          h1_pid = list1[a].pid;
          h1_mass = m1;
          h2_pid = sign2 * list2[b].pid;
          h2_mass = m2;
          return total;
        }
      }
    }
  }

  return total;
}

// ---------------------------------------------------------------------------
// Main Decay Function

__global__ void decay_clusters(event* events, cluster_list* cls,
                               hadronisation* had, int n) {
  /**
   * @brief Decay each cluster into two hadrons, or one if it reaches no pair.
   *
   * Baryonic with probability pwt[3] / (1 + pwt[3]), falling back to mesons;
   * the pair is drawn with probability W / sum(W) (see choose_decay).
   *
   * @param events the events to process
   * @param cls the per-event cluster lists
   * @param had the hadronisation object
   * @param n the number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  cluster_list& cl = cls[idx];
  // Check for overflow
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  // In case clusters cannot decay into two hadrons, set their mass to
  // the mass of the lightest meson (q1,q2) forms
  fallback_reshuffling(ev, cl);

  // loop over all clusters
  int n_cl = cl.get_size();

  for (int i = 0; i < n_cl; i++) {
    cluster c = cl.get_cluster(i);

    // constituents
    int i1 = c.get_i1();
    int i2 = c.get_i2();
    int p1_pid = ev.get_particle(i1).get_pid();  // quark  (positive)
    int p2_pid = ev.get_particle(i2).get_pid();  // antiquark (negative)

    vec4 cl_mom = c.get_mom();
    double M = cl_mom.m();

    // Category: baryonic with probability pwt[3] / (1 + pwt[3]), and only
    // where the cluster reaches a baryon pair at all. The draw is made either
    // way, so the random stream is the same as it would be without the question
    double rho_cat = ev.gen_random();
    bool baryonic = (rho_cat * (1.0 + had->pwt[3]) < had->pwt[3]) &&
                    can_produce_baryon_pair(p1_pid, p2_pid, M);

    // Whether the category opens anything at all. A baryonic draw has already
    // asked; a mesonic one reaches a pair exactly when the cluster is above
    // the lightest one, so this is the walk's own answer, arrived at without
    // walking. A cluster that opens nothing takes the single hadron below
    bool reachable = baryonic || M > lightest_meson_pair_mass(p1_pid, p2_pid);

    if (reachable) {
      // The pair drawn, and the masses it decays with
      int h1_pid = 0;
      int h2_pid = 0;
      double h1_mass = 0.;
      double h2_mass = 0.;

      // Total the weight of every pair this category opens
      double total = choose_decay(had->pwt, p1_pid, p2_pid, M, baryonic, -1.,
                                  h1_pid, h1_mass, h2_pid, h2_mass);

      // Draw one pair with probability W / sum(W)
      double rho_pair = ev.gen_random();
      choose_decay(had->pwt, p1_pid, p2_pid, M, baryonic, rho_pair * total,
                   h1_pid, h1_mass, h2_pid, h2_mass);

      // A category that opens a pair always selects one: the draw is below the
      // total, and the two walks accumulate the same pairs in the same order.
      // A miss means those walks disagreed, which no correct run reaches
      if (h1_pid == 0) {
        printf("decay_clusters: no pair selected for a reachable cluster\n");
      } else {
        // Decay the cluster into the two hadrons
        vec4 h1, h2;
        double rho_1 = ev.gen_random();
        double rho_2 = ev.gen_random();
        one_to_two_decay(cl_mom, h1_mass, h2_mass, h1, h2, rho_1, rho_2);

        ev.set_particle(i1, particle(h1_pid, h1, 0, 0));
        ev.set_particle(i2, particle(h2_pid, h2, 0, 0));

        // Continue
        continue;
      }
    }

    // Single hadron fallback (q1, q2) -> hadron
    candidate hadron = lightest_meson(p1_pid, p2_pid);
    ev.set_particle(i1, particle(hadron.pid, cl_mom, 0, 0));
    ev.set_particle(i2, particle(0, vec4(), 0, 0));
  }

  // Drop the slots the single-hadron fallback emptied. The cluster list is
  // spent by now, so nothing holds an index across this.
  ev.compact();
}
