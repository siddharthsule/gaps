#include "shower.cuh"

// -----------------------------------------------------------------------------
// constructor

__device__ void shower::setup(double t_c, double as_max) {
  /**
   * @brief construct the shower; the pdf and alpha_s objects stay on the host
   *
   * @param t_c the shower cutoff scale
   * @param as_max the maximum value of alpha_s, the veto overestimate
   */
  this->t_c = t_c;
  this->as_max = as_max;
}

// kernel to set up the shower object on the device
__global__ void shower_setup_kernel(shower* sh, double t_c, double as_max) {
  /**
   * @brief Set up the shower object on the device
   *
   * @param sh The shower object
   * @param t_c The cutoff scale
   * @param as_max The maximum value of alpha_s
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= 1) return;
  // ---------------------------------------------
  sh->setup(t_c, as_max);
}

// -----------------------------------------------------------------------------
// preparing the shower

__global__ void prep_shower(event* events, bool nlo_matching, int n) {
  /**
   * @brief Prepares the shower for the event
   *
   * @param events The events to prepare
   * @param nlo_matching Whether NLO matching provides the first emission
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  event& ev = events[idx];
  // Nothing to shower into if the record is already full
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  // NLO Matching does the first emission and sets the shower scale
  // to the first emission pT. If NLO Matching is off, we find the
  // smallest pT in the event and set the shower scale to that.
  if (!nlo_matching) {
    // set the starting shower scale
    double t_start = 10000000.;
    for (int i = 0; i < ev.get_size(); i++) {
      for (int j = i + 1; j < ev.get_size(); j++) {
        // Skip non-partons
        if (!ev.get_particle(i).is_parton() ||
            !ev.get_particle(j).is_parton()) {
          continue;
        }

        // Check if ij and k are colour connected
        if (!ev.get_particle(i).is_color_connected(ev.get_particle(j))) {
          continue;
        }

        // Get the invariant mass squared
        double t =
            (ev.get_particle(i).get_mom() + ev.get_particle(j).get_mom()).m2();

        // Check if minimum
        if (t < t_start) {
          t_start = t;
        }
      }
    }

    // Count colour connections: gluons have 2, quarks/antiquarks have 1
    double c_start = 0;
    for (int i = 0; i < ev.get_size(); i++) {
      if (!ev.get_particle(i).is_parton()) {
        continue;
      }
      // Gluons (pid=21) have both colour and anticolour
      if (abs(ev.get_particle(i).get_pid()) == 21) {
        c_start += 2.;
      } else {
        // Quarks and antiquarks have one colour connection
        c_start += 1.;
      }
    }
    // Each colour line connects two partons, so divide by 2
    c_start /= 2.;

    ev.set_shower_t(t_start);
    ev.set_shower_c(c_start);
  }
}

// -----------------------------------------------------------------------------

__global__ void select_winner_split_func(shower* shower, event* events,
                                         int* active_idx, int n,
                                         double* winner) {
  /**
   * @brief Select the winner (highest transverse momentum) trial splitting
   *
   * This kernel takes about half of the shower time.
   *
   * @param shower The shower object
   * @param events The events to run the shower on
   * @param active_idx The event index held by each active slot
   * @param n The number of active slots
   * @param winner The array to store the winner emission data
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[active_idx[idx]].has_shower_ended()) return;
  event& ev = events[active_idx[idx]];
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  // default values
  double win_tt = shower->t_c;
  int win_sf = 0;  // 0 = no splitting
  int win_ij = 0;
  int win_k = 0;
  double win_sijk = 0.;
  double win_zm = 0.;
  double win_zp = 0.;

  for (int ij = 0; ij < ev.get_size(); ij++) {
    /**
     * This double for loop is quite expensive, but thanks to QCD, we know that
     * each parton can only have up to two colour connected partners, so we can
     * break the inner loop after we have checked these partners.
     */
    int partners_checked = 0;

    for (int k = 0; k < ev.get_size(); k++) {
      // sanity check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // Skip non-partons
      if (!ev.get_particle(ij).is_parton() || !ev.get_particle(k).is_parton()) {
        continue;
      }

      // break this loop if we have checked:
      // For a quark - only 1 partner
      // For a gluon - 2 partners
      if (ev.get_particle(ij).get_pid() == 21 && partners_checked == 2) {
        break;
      } else if (ev.get_particle(ij).get_pid() != 21 && partners_checked == 1) {
        break;
      }

      // need to check if ij and k are colour connected
      if (!ev.get_particle(ij).is_color_connected(ev.get_particle(k))) {
        continue;
      }

      // Increment partners_checked
      partners_checked++;

      // get the invariant mass squared of the dipole
      double sijk =
          (ev.get_particle(ij).get_mom() + ev.get_particle(k).get_mom()).m2();

      // get the splitting functions for the current partons
      // sf_codes is an array of possible splitting functions
      int sf_codes[11];
      shower->generate_possible_splittings(
          ev.get_particle(ij).get_pid(), ev.get_particle(k).get_pid(),
          ev.get_particle(ij).is_initial(), ev.get_particle(k).is_initial(),
          sf_codes);

      // codes instead of object oriented approach! See splittings.cu
      for (int sf : sf_codes) {
        // When a null code is encountered, we have reached the end of the
        // possible splittings, we can break out of the loop
        if (sf == -1) {
          break;
        }

        // check if either parton has eta < eta_min (usually 1e-5)
        if (shower->is_ii(sf) && (ev.get_particle(ij).get_eta() < 1e-5 ||
                                  ev.get_particle(k).get_eta() < 1e-5)) {
          continue;
        }

        // phase space limits
        double zm, zp;

        // If FI, send eta_k for the limits
        // If IF/II, send eta_ij
        double eta = shower->is_fi(sf) ? ev.get_particle(k).get_eta()
                                       : ev.get_particle(ij).get_eta();
        shower->get_boundaries(zm, zp, sijk, eta, sf);

        if (zm < 0. || zp > 1. || zm > zp) {
          continue;
        }

        // Calculate the integrated overestimate
        double pdf_max = shower->get_pdf_max(sf, ev.get_particle(ij).get_eta());
        double j0_max = shower->is_ff(sf) ? 1. : 2.;
        double c = shower->as_max / (2. * M_PI) *
                   shower->sf_integral(zm, zp, sf) * j0_max * pdf_max;

        // calculate the evolution variable
        double tt;

        // Final State Radiation - As always
        if (shower->is_fsr(sf)) {
          // t = T * random^(1/c)
          tt = ev.get_shower_t() * pow(ev.gen_random(), 1. / c);
        }

        // Initial State Radiation - Need to account for quark masses for PDF
        // Ranges (Avoid evolving to below quark mass)
        else {
          // Get the quark mass squared
          int fl = abs(ev.get_particle(ij).get_pid());
          double mq2 = fl == 5 ? mb * mb : (fl == 4 ? mc * mc : 0.);

          // Evolve from (t - m2) not t for this case
          tt = (ev.get_shower_t() - mq2) * pow(ev.gen_random(), 1. / c) + mq2;

          // Check if tt <= mq2: Need some numerical error protection
          if (tt - mq2 < 1e-9) {
            continue;
          }
        }

        // If g->bb or g->cc, check if tt is above the quark mass threshold, so
        // the hadronisation reshuffler need not push c/b onto their masses
        if (shower->is_g2qqbar(sf) || shower->is_g2qbarq(sf)) {
          if ((shower->get_splitting_flavour(sf) == 5 && tt < mb2) ||
              (shower->get_splitting_flavour(sf) == 4 && tt < mc2)) {
            continue;
          }
        }

        // For FI, IF and II, skip events where q2 (=tt) is less than pdf limit
        if (!shower->is_ff(sf) && tt < (pdf_q_min * pdf_q_min)) {
          continue;
        }

        // check if tt is greater than the current winner
        if (tt > win_tt) {
          win_tt = tt;
          win_sf = sf;
          win_ij = ij;
          win_k = k;
          win_sijk = sijk;
          win_zm = zm;
          win_zp = zp;
        }
      }
    }
  }

  // set the new shower t
  ev.set_shower_t(win_tt);

  // Also generate z, y and phi
  double z = shower->sf_generate_z(win_zm, win_zp, ev.gen_random(), win_sf);
  double y = shower->calculate_y(win_tt, z, win_sijk, win_sf);
  double phi = 2. * M_PI * ev.gen_random();

  // IF: z, pt -> x, u (here, z, y)
  if (shower->is_if(win_sf)) {
    double z0 = z;
    double ratio = win_tt / win_sijk;
    double frac2 =
        4 * ratio * z0 * (1. - z0) / ((1. - z0 + ratio) * (1. - z0 + ratio));
    z = (1. - z0 + ratio) / (2 * ratio) * (1. - sqrt(1. - frac2));
    y = (z * ratio) / (1. - z0);
  }

  // II: z, pt -> x, v (here, z, y)
  else if (shower->is_ii(win_sf)) {
    double z0 = z;
    double ratio = win_tt / win_sijk;
    z = z0 * (1. - z0) / (1. - z0 + ratio);
    y = (z * ratio) / (1. - z0);
  }

  // Set the winner variables (sf, ij, k, sijk, z, y, phi)
  winner[7 * idx] = static_cast<double>(win_sf);
  winner[7 * idx + 1] = static_cast<double>(win_ij);
  winner[7 * idx + 2] = static_cast<double>(win_k);
  winner[7 * idx + 3] = win_sijk;
  winner[7 * idx + 4] = z;
  winner[7 * idx + 5] = y;
  winner[7 * idx + 6] = phi;
}

// -----------------------------------------------------------------------------

__global__ void check_cutoff(event* events, int* active_idx, shower* shower,
                             int* d_completed, int n) {
  /**
   * @brief Check if the shower has ended
   *
   * @param events The events to run the shower on
   * @param active_idx The event index held by each active slot
   * @param shower The shower object
   * @param d_completed The number of completed events
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[active_idx[idx]].has_shower_ended()) return;
  event& ev = events[active_idx[idx]];
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  // end shower if t < cutoff, written !(t > t_c) as in the literature
  if (!(ev.get_shower_t() > shower->t_c)) {
    ev.shower_has_ended(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events

    return;
  }
}

// -----------------------------------------------------------------------------

// PDF ratio: xf is evaluated for ij and i (see pdf.cuh), and the ratio is
// taken in the veto_alg kernel

// -----------------------------------------------------------------------------

__global__ void veto_alg(shower* shower, alpha_s* as, event* events,
                         int* active_idx, int n, double* xf_a, double* xf_b,
                         bool* accept_emission, double* winner) {
  /**
   * @brief The veto algorithm for the shower
   *
   * @param shower The shower object
   * @param as The alpha_s object
   * @param events The events to run the shower on
   * @param active_idx The event index held by each active slot
   * @param n The number of active slots
   * @param xf_a The PDF of the parton after emission
   * @param xf_b The PDF of the parton before emissions
   * @param accept_emission The array to store the acceptance of the emission
   * @param winner The array to store the winner emission data
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[active_idx[idx]].has_shower_ended()) return;
  event& ev = events[active_idx[idx]];
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  // set to false, only set to true if accepted
  accept_emission[idx] = false;

  // Get the shower evolution variable
  double t = ev.get_shower_t();

  // Get the winner variables (sf, ij, j, sijk, z, y, phi)
  int sf = static_cast<int>(winner[7 * idx]);
  int ij = static_cast<int>(winner[7 * idx + 1]);
  int k = static_cast<int>(winner[7 * idx + 2]);
  double sijk = winner[7 * idx + 3];
  double z = winner[7 * idx + 4];
  double y = winner[7 * idx + 5];

  // calculate z (as in the sampled variable)
  double z0;
  if (shower->is_isr(sf)) {
    z0 = 1. - (z * (t / sijk) / y);
  } else {
    // FF and FI: can use z directly
    z0 = z;
  }

  // Check for Phase Space
  if (!(shower->check_phase_space(z, y, sf))) {
    return;
  }

  // Get PDF Ratio and PDF Max for FI and IF/II
  double pdf_ratio(1.), pdf_max(1.);

  if (!shower->is_ff(sf)) {
    // Check Momentum Fraction
    if (!(shower->check_mom_frac(sf, ev.get_particle(ij).get_pid(),
                                 ev.get_particle(k).get_pid(),
                                 ev.get_particle(ij).get_eta(),
                                 ev.get_particle(k).get_eta(), z, y))) {
      return;
    }

    // Ensure PDFs are valid
    if (isnan(xf_a[idx]) || isnan(xf_b[idx]) || isinf(xf_a[idx]) ||
        isinf(xf_b[idx]) || xf_a[idx] <= 0. || xf_b[idx] <= 0.) {
      return;
    }

    // Calculate the ratio and handle division by zero
    pdf_ratio = xf_a[idx] / xf_b[idx];

    // cancel emission if pdf_ratio is nan, inf, or less than 0
    if (isnan(pdf_ratio) || isinf(pdf_ratio) || (pdf_ratio <= 0.)) {
      return;
    }

    // If PDF Ratio is too large, veto
    if (pdf_ratio > 1e6) {
      return;
    }

    // LHAPDF gives xf(x, q2), so the ratio is 1/z * f(x/z, q2) / f(x, q2);
    // the factor of z removes the 1/z here rather than in the jacobian
    pdf_ratio *= z;

    // Mutliply by (t - m2) / t for ISR to account for quark masses
    if (shower->is_isr(sf)) {
      int fl = abs(ev.get_particle(ij).get_pid());
      pdf_ratio *= (t - (fl == 5 ? mb * mb : (fl == 4 ? mc * mc : 0.))) / t;
    }

    // Get PDF Max
    pdf_max = shower->get_pdf_max(sf, ev.get_particle(ij).get_eta());
  }

  // Jacobian
  double jacobian = shower->get_jacobian(z, y, sf) * pdf_ratio;
  double j0_max = shower->is_ff(sf) ? 1. : 2.;
  double jmaxtot = j0_max * pdf_max;

  // Splitting Function Value and Estimate
  double value = shower->sf_value(z, y, sf);
  double estimate = shower->sf_estimate(z0, sf);

  // veto algorithm
  double f = (*as)(t) / (2. * M_PI) * value * jacobian;
  double g = shower->as_max / (2. * M_PI) * estimate * jmaxtot;

  // Check for Negative f
  if (f < 0.) {
    return;
  }

  // Check for f > g
  if (f > g) {
    return;
  }

  // Accept / Veto
  if (ev.gen_random() < f / g) {
    accept_emission[idx] = true;
  }
}

// -----------------------------------------------------------------------------

// do splitting
__global__ void do_splitting(shower* shower, event* events, int* active_idx,
                             int n, bool* accept_emission, double* winner) {
  /**
   * @brief Do the splitting for the shower
   *
   * @param shower The shower object
   * @param events The events to run the shower on
   * @param active_idx The event index held by each active slot
   * @param n The number of active slots
   * @param accept_emission The array to store the acceptance of the emission
   * @param winner The array to store the winner emission data
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[active_idx[idx]].has_shower_ended()) return;
  event& ev = events[active_idx[idx]];
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  // Do not run if the emission was not accepted by the veto algorithm
  if (!accept_emission[idx]) {
    return;
  }

  // Get the shower evolution variable
  double t = ev.get_shower_t();

  // Get the winner variables (sf, ij, j, sijk, z, y, phi)
  int sf = static_cast<int>(winner[7 * idx]);
  int ij = static_cast<int>(winner[7 * idx + 1]);
  int k = static_cast<int>(winner[7 * idx + 2]);
  double z = winner[7 * idx + 4];
  double y = winner[7 * idx + 5];
  double phi = winner[7 * idx + 6];

  // get the flavours
  int flavs[3];
  shower->sf_to_flavs(sf, flavs);

  // pi, pj, pk, pijt, pkt and kt
  vec4 moms[6] = {vec4(), vec4(), vec4(), vec4(), vec4(), vec4()};
  shower->make_kinematics(moms, z, y, phi, ev.get_particle(ij).get_mom(),
                          ev.get_particle(k).get_mom(), sf);

  // calculate the colours
  int colij[2] = {ev.get_particle(ij).get_col(),
                  ev.get_particle(ij).get_acol()};
  int colk[2] = {ev.get_particle(k).get_col(), ev.get_particle(k).get_acol()};
  int coli[2] = {0, 0};
  int colj[2] = {0, 0};
  shower->make_colours(ev.get_shower_c(), sf, flavs, colij, colk, coli, colj,
                       ev.gen_random());

  // modify splitter
  ev.set_particle_pid(ij, flavs[1]);
  ev.set_particle_mom(ij, moms[0]);
  ev.set_particle_col(ij, coli[0]);
  ev.set_particle_acol(ij, coli[1]);
  if (shower->is_isr(sf)) {
    ev.set_particle_eta(ij, ev.get_particle(ij).get_eta() / z);
  }

  // modify recoiled spectator
  ev.set_particle_mom(k, moms[2]);
  if (shower->is_fi(sf)) {
    ev.set_particle_eta(k, ev.get_particle(k).get_eta() / y);
  }

  // add emitted parton (return = overflowed)
  particle em = particle(flavs[2], moms[1], colj[0], colj[1]);
  if (!ev.add_emission(em)) return;

  // II Only - Lorentz Boost the new final state
  if (shower->is_ii(sf)) {
    shower->ii_boost_after_emission(ev, moms);
  }

  return;
}

// -----------------------------------------------------------------------------

__global__ void check_too_many_particles(event* events, int* active_idx,
                                         int n_emissions_max,
                                         int* d_completed, int n) {
  /**
   * @brief Check if the event has too many particles
   *
   * Ends overflowed events rather than skipping them, as this kernel owns
   * their d_completed count, which the shower loop waits on.
   *
   * @param events The events to run the shower on
   * @param active_idx The event index held by each active slot
   * @param n_emissions_max The maximum number of emissions
   * @param d_completed The number of completed events
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[active_idx[idx]].has_shower_ended()) return;
  event& ev = events[active_idx[idx]];
  // ---------------------------------------------

  // An overflowed event has nowhere left to emit, whatever its scale
  if (ev.get_overflowed()) {
    ev.shower_has_ended(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events
    return;
  }

  // limit to max particles
  if (ev.get_size() == min(max_particles, ev.get_hard() + n_emissions_max)) {
    // The record running out is an overflow, and the event is dropped
    if (ev.get_size() >= max_particles) ev.set_overflowed();
    ev.shower_has_ended(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events
    return;
  }
}

// -----------------------------------------------------------------------------

struct is_active_event {
  /**
   * @brief Function object to check if an event is still showering
   *
   * @param a_idx The index of the event in the event array
   * @return true if the shower has not ended
   */
  event* events;

  __device__ bool operator()(int a_idx) const {
    return !events[a_idx].has_shower_ended();
  }
};

// -----------------------------------------------------------------------------

void run_shower(thrust::device_vector<event>& dv_events, const params& p,
                int blocks) {
  /**
   * @brief Run the shower on the events
   *
   * @param dv_events The events to run the shower on
   * @param p The run parameters (cutoff, alpha_s, pdf, nlo, partitioning)
   * @param blocks The number of thread blocks to launch the kernels with
   */

  // number of events - can get from d_events.size()
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // set up the device alpha_s calculator
  alpha_s* d_as;
  cudaMalloc(&d_as, sizeof(alpha_s));
  as_setup_kernel<<<1, 1>>>(d_as, p.asmz, (p.fixed_as ? 0 : 2), p.use_cmw);
  sync_gpu_and_check("as_setup_kernel");

  // Calculate as_max = as(t_c)
  double* d_as_max;
  cudaMalloc(&d_as_max, sizeof(double));
  as_value<<<1, 1>>>(d_as, d_as_max, p.t_c);
  sync_gpu_and_check("as_value");
  double as_max;
  cudaMemcpy(&as_max, d_as_max, sizeof(double), cudaMemcpyDeviceToHost);

  // set up the shower
  shower* d_shower;
  cudaMalloc(&d_shower, sizeof(shower));
  shower_setup_kernel<<<1, 1>>>(d_shower, p.t_c, as_max);
  sync_gpu_and_check("shower_setup_kernel");

  // set up the pdf evaluator
  pdf_wrapper pdf(p.showerpdf);

  // Winner variables (sf, ij, k, sijk, z, y, phi), 7 per event, all stored
  // as doubles (static_cast<int> for sf, ij, k); t and c live in the event
  thrust::device_vector<double> dv_winner(7 * n_events, 0.0);
  double* d_winner = thrust::raw_pointer_cast(dv_winner.data());

  // pdf ratio
  thrust::device_vector<double> dv_x_a(n_events, 0.0);
  thrust::device_vector<double> dv_x_b(n_events, 0.0);
  thrust::device_vector<double> dv_q2(n_events, 0.0);
  thrust::device_vector<double> dv_xf_a(n_events, 0.0);
  thrust::device_vector<double> dv_xf_b(n_events, 0.0);
  thrust::device_vector<double> dv_ratio(n_events, 0.0);
  thrust::device_vector<int> dv_flavours_a(n_events, 0);
  thrust::device_vector<int> dv_flavours_b(n_events, 0);
  double* d_x_a = thrust::raw_pointer_cast(dv_x_a.data());
  double* d_x_b = thrust::raw_pointer_cast(dv_x_b.data());
  double* d_q2 = thrust::raw_pointer_cast(dv_q2.data());
  double* d_xf_a = thrust::raw_pointer_cast(dv_xf_a.data());
  double* d_xf_b = thrust::raw_pointer_cast(dv_xf_b.data());
  double* d_ratio = thrust::raw_pointer_cast(dv_ratio.data());
  int* d_flavours_a = thrust::raw_pointer_cast(dv_flavours_a.data());
  int* d_flavours_b = thrust::raw_pointer_cast(dv_flavours_b.data());

  // veto outcome
  thrust::device_vector<bool> dv_accept_emission(n_events, false);
  bool* d_accept_emission = thrust::raw_pointer_cast(dv_accept_emission.data());

  /**
   * Active event indices
   * --------------------
   * dv_active_idx holds the positions in dv_events of the events still being
   * showered. Kernels run over slots [0, n), reach their event through
   * active_idx[idx], and index the per-event buffers above by slot. Each cycle
   * the list is compacted, so finished events are no longer visited and the
   * event records never move.
   */
  thrust::device_vector<int> dv_active_idx(n_events);
  thrust::device_vector<int> dv_next_idx(n_events);
  thrust::sequence(dv_active_idx.begin(), dv_active_idx.end());
  int* d_active_idx = thrust::raw_pointer_cast(dv_active_idx.data());

  // ---------------------------------------------------
  // Analysis Variables

  // allocate device memory for completed events counter
  int* d_completed;
  cudaMalloc(&d_completed, sizeof(int));
  cudaMemset(d_completed, 0, sizeof(int));

  // ---------------------------------------------------------------------------
  // prepare the shower

  debug_msg("running @prep_shower");
  prep_shower<<<blocks, p.threads>>>(d_events, p.nlo, n_events);
  sync_gpu_and_check("prep_shower");

  // ---------------------------------------------------------------------------
  // run the shower

  // number of completed events and cycles
  int completed = 0;
  int cycle = 0;

  // (Varying) kernel size
  int n = n_events;

  // An event whose starting scale is already at or below the cutoff has no
  // shower to run. End it here, before select_winner draws for it: the CPU
  // loop never enters for such an event, so its PRNG has to be left where
  // the CPU leaves it for the hadronisation to follow the same stream. This
  // occurs with NLO matching
  debug_msg("running @check_cutoff (initial)");
  check_cutoff<<<blocks, p.threads>>>(d_events, d_active_idx, d_shower,
                                      d_completed, n);
  sync_gpu_and_check("check_cutoff (initial)");
  cudaMemcpy(&completed, d_completed, sizeof(int), cudaMemcpyDeviceToHost);

  while (completed < n_events) {
    // run all the kernels here...

    // -------------------------------------------------------------------------
    // check if there are too many particles (do first in case of H event)

    debug_msg("running @check_too_many_particles");
    check_too_many_particles<<<blocks, p.threads>>>(
        d_events, d_active_idx, p.n_emissions_max, d_completed, n);
    sync_gpu_and_check("check_too_many_particles");

    // -------------------------------------------------------------------------
    // select the winner kernel

    debug_msg("running @select_winner_split_func");
    select_winner_split_func<<<blocks, p.threads>>>(d_shower, d_events,
                                                    d_active_idx, n, d_winner);
    sync_gpu_and_check("select_winner_split_func");

    // -------------------------------------------------------------------------
    // check cutoff

    debug_msg("running @check_cutoff");
    check_cutoff<<<blocks, p.threads>>>(d_events, d_active_idx, d_shower,
                                        d_completed, n);
    sync_gpu_and_check("check_cutoff");

    // -------------------------------------------------------------------------
    // calculate pdf of ij and i using LHAPDF

    // Skip for LEP
    if (p.process != 1) {
      debug_msg("running @setup_pdfratio");
      setup_pdfratio<<<blocks, p.threads>>>(d_shower, d_events, d_active_idx, n,
                                            d_flavours_a, d_flavours_b, d_x_a,
                                            d_x_b, d_q2, d_winner);
      sync_gpu_and_check("setup_pdfratio");

      pdf.evaluate(d_flavours_a, d_x_a, d_q2, d_xf_a, n, blocks, p.threads);
      pdf.evaluate(d_flavours_b, d_x_b, d_q2, d_xf_b, n, blocks, p.threads);
    }

    // -------------------------------------------------------------------------
    // veto algorithm

    debug_msg("running @veto_alg");
    veto_alg<<<blocks, p.threads>>>(d_shower, d_as, d_events, d_active_idx, n,
                                    d_xf_a, d_xf_b, d_accept_emission,
                                    d_winner);
    sync_gpu_and_check("veto_alg");

    // -------------------------------------------------------------------------
    // splitting algorithm

    debug_msg("running @do_splitting");
    do_splitting<<<blocks, p.threads>>>(d_shower, d_events, d_active_idx, n,
                                        d_accept_emission, d_winner);
    sync_gpu_and_check("do_splitting");

    // -------------------------------------------------------------------------
    // import the number of completed events

    cudaMemcpy(&completed, d_completed, sizeof(int), cudaMemcpyDeviceToHost);
    cycle++;

    // -------------------------------------------------------------------------

    // Update the active indices
    // Stop at 25k, below this, the overhead of partitioning outweighs the
    // benefits, increasing execution time.
    if (p.do_partitioning && n > 25000) {
      // Copy the indices of the still-active events in the first n slots
      // into dv_next_idx, keeping their order
      auto next_end =
          thrust::copy_if(dv_active_idx.begin(), dv_active_idx.begin() + n,
                          dv_next_idx.begin(), is_active_event{d_events});

      // Update n to reflect the number of incomplete events
      n = static_cast<int>(next_end - dv_next_idx.begin());

      // The compacted indices become the active set
      dv_active_idx.swap(dv_next_idx);
      d_active_idx = thrust::raw_pointer_cast(dv_active_idx.data());
    }

    // -------------------------------------------------------------------------
    // show progress

    std::cerr << "\rCompleted Events: " << completed << "/" << n_events
              << std::flush;

    // -------------------------------------------------------------------------
  }
  std::cout << std::endl;

  // ---------------------------------------------------------------------------

  // free the memory
  cudaFree(d_shower);
  cudaFree(d_as);
  cudaFree(d_completed);
}