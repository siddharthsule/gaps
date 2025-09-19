#include "shower.cuh"

// -----------------------------------------------------------------------------
// Kernel to setup the arrays for pdf ratio calculation

__global__ void setup_pdfratio(shower *shower, event *events, int n,
                               int *flavours_a, int *flavours_b, double *x_a,
                               double *x_b, double *q2, double *winner) {
  /**
   * @brief Setup the arrays for the pdf ratio calculation
   *
   * @param shower: Shower object
   * @param events: Array of event records
   * @param n: Number of events
   * @param flavours_a: Array of flavours for parton a (after)
   * @param flavours_b: Array of flavours for parton b (before)
   * @param x_a: Array of momentum fractions for parton a
   * @param x_b: Array of momentum fractions for parton b
   * @param q2: Array of evolution variables
   * @param winner: Array of winner variables
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[idx].has_shower_ended()) return;
  event &ev = events[idx];
  // ---------------------------------------------

  // Get the shower evolution variable
  double t = ev.get_shower_t();

  // Get the winner variables (sf, ij, j, sijk, z, y, phi)
  int sf = static_cast<int>(winner[7 * idx]);
  int ij = static_cast<int>(winner[7 * idx + 1]);
  int k = static_cast<int>(winner[7 * idx + 2]);
  double z = winner[7 * idx + 4];

  // No need to fill if FF
  if (shower->is_ff(sf)) return;

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = shower->get_splitting_case(sf);

  // 0 = q->qg, 1 = q->gq, 2 = g->gg, 3 = g->qq
  int splitting_type = shower->get_splitting_type(sf);

  // 0 = particle, 1 = antiparticle
  int is_em_antip = shower->is_emitter_antiparticle(sf);

  // 1 = d, 2 = u, 3 = s, 4 = c, 5 = b; 0 = gluon
  int splitting_flavour = shower->get_splitting_flavour(sf);

  // Get the parton data
  int ij_pid = ev.get_particle(ij).get_pid();
  int k_pid = ev.get_particle(k).get_pid();
  double ij_eta = ev.get_particle(ij).get_eta();
  double k_eta = ev.get_particle(k).get_eta();

  // NOTE: a = after, b = before

  switch (splitting_case) {
    case 1:
      flavours_a[idx] = k_pid;
      flavours_b[idx] = k_pid;
      x_a[idx] = k_eta / z;
      x_b[idx] = k_eta;
      q2[idx] = t;
      break;

    case 2:
    case 3:

      switch (splitting_type) {
        // q -> qg
        case 0:
          flavours_a[idx] = ij_pid;
          flavours_b[idx] = ij_pid;
          break;

        // q -> gq
        case 1:
          flavours_a[idx] = 21;
          flavours_b[idx] = ij_pid;
          break;

        // g -> gg
        case 2:
          flavours_a[idx] = ij_pid;
          flavours_b[idx] = ij_pid;
          break;

        // g -> q qbar (backwards is qbar to g qbar)
        case 3:
          flavours_a[idx] = -1 * splitting_flavour;
          flavours_b[idx] = 21;
          break;

        // g -> qbar q (backwards is q to g q)
        case 4:
          flavours_a[idx] = splitting_flavour;
          flavours_b[idx] = 21;
          break;

        default:
          break;
      }

      x_a[idx] = ij_eta / z;
      x_b[idx] = ij_eta;
      q2[idx] = t;
      break;

    default:
      break;
  }
  return;
}

// -----------------------------------------------------------------------------
// Check if the momentum fraction is valid before calculating

__device__ bool shower::check_mom_frac(int sf, int ij_pid, int k_pid,
                                       double ij_eta, double k_eta,
                                       double z) const {
  /**
   * @brief Check if the momentum fraction is valid before calculating the pdf
   * ratio
   *
   * @param sf Splitting Function Code
   * @param ij_pid Particle ID of the emitter
   * @param k_pid Particle ID of the spectator
   * @param ij_eta Momentum Fraction of the emitter
   * @param k_eta Momentum Fraction of the spectator
   * @param z Splitting Variable
   * @return true: If the momentum fraction is valid
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  switch (splitting_case) {
    // FI splittings: z and x (here, y)
    case 1:
      return (k_eta < z) && (k_eta > pdf_x_min) && (k_eta < pdf_x_max) &&
             (k_eta / z > pdf_x_min) && (k_eta / z < pdf_x_max);
      break;

    // IF splittings: x (here, z) and u (here, y)
    case 2:
    case 3:
      return (ij_eta < z) && (ij_eta > pdf_x_min) && (ij_eta < pdf_x_max) &&
             (ij_eta / z > pdf_x_min) && (ij_eta / z < pdf_x_max);
      break;
  }
  return false;
}

// -----------------------------------------------------------------------------
// Get the appropriate pdf_max value

__device__ double shower::get_pdf_max(int sf, double ij_eta) const {
  /**
   * @brief Get the maximum pdf ratio for the given splitting
   *
   * @param sf Splitting Function Code
   * @param ij_eta Momentum Fraction of the emitter
   * @return double: PDF Max
   */

  // No pdf_max for FF splittings
  if (get_splitting_case(sf) == 0) {
    return 1.;
  }

  // IF/II q -> gq type splittings
  if ((get_splitting_case(sf) > 1) && (get_splitting_type(sf) == 1)) {
    // qbar -> g qbar
    if (is_emitter_antiparticle(sf)) {
      switch (get_splitting_flavour(sf)) {
        case 5:
          return fmax(150., 700. * ij_eta);
        case 4:
          return fmax(100., 450. * ij_eta);
        case 3:
          return fmax(50., 600. * ij_eta);
        default:
          return 150.;
      }
    }

    // q -> g q
    else {
      switch (get_splitting_flavour(sf)) {
        case 5:
          return fmax(150., 650. * ij_eta);
        case 4:
          return fmax(100., 450. * ij_eta);
        case 3:
          return fmax(50., 600. * ij_eta);
        default:
          return 150.;
      }
    }
  }

  // Iniital State g -> u ubar or g -> ubar u
  // Small unless eta > 0.1 (gluon wants to become an up quark for proton!)
  else if ((sf == 3302) || (sf == 3402) || (sf == 2302) || (sf == 2402)) {
    if (ij_eta > 0.1) {
      return 10.;
    }
    return 2.;
  }

  // All other Splittings
  else {
    return 2.;
  }
}