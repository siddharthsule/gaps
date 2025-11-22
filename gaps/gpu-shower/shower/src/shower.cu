#include "shower.cuh"

// -----------------------------------------------------------------------------
// constructor

__device__ void shower::setup(double root_s, double t_c, double as_max) {
  this->e_proton = root_s / 2.;
  this->t_c = t_c;
  this->as_max = as_max;
  this->j0_max = 2.;
}

// kernel to set up the matrix object on the device
__global__ void shower_setup_kernel(shower* sh, double root_s, double t_c,
                                    double as_max) {
  /**
   * @brief Set up the shower object on the device
   *
   * @param sh The shower object
   * @param as The alpha_s object
   * @param root_s The root s energy
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= 1) return;
  // ---------------------------------------------
  sh->setup(root_s, t_c, as_max);
}

// -----------------------------------------------------------------------------
// preparing the shower

__global__ void prep_shower(event* events, bool nlo_matching, int n) {
  /**
   * @brief Prepares the shower for the event
   *
   * @param events The events to prepare
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  event& ev = events[idx];
  // ---------------------------------------------

  // NLO Matching does the first emission and sets the shower scale
  // to the first emission pT. If NLO Matching is off, we find the
  // smallest pT in the event and set the shower scale to that.
  if (!nlo_matching) {
    // set the starting shower scale
    double t_start = 10000000.;
    for (int i = 0; i < ev.get_size(); i++) {
      for (int j = 0; j < ev.get_size(); j++) {
        if (i == j) {
          continue;
        }

        if (!ev.get_particle(i).is_parton() ||
            !ev.get_particle(j).is_parton()) {
          continue;
        }

        double t =
            (ev.get_particle(i).get_mom() + ev.get_particle(j).get_mom()).m2();

        if (t < t_start) {
          t_start = t;
        }
      }
    }

    ev.set_shower_t(t_start);
    ev.set_shower_c(1);
  }
}

// -----------------------------------------------------------------------------

__global__ void select_winner_split_func(shower* shower, event* events, int n,
                                         double* winner) {
  /**
   * @brief Select the winner splitting in the event
   *
   * This function generates the highest transverse momentum splitting for every
   * possible dipoles in the event. It then chooses the winner emission from the
   * generated splittings, by picking the one with the highest transverse
   * momentum. This winner emission is then used in the veto step
   *
   * when you profile the code, you will notice that this is the process that
   * takes up half of the shower time. this method below is a first attempt at
   * parallelizing the process.
   *
   * @param shower The shower object
   * @param events The events to run the shower on
   * @param n The number of events
   * @param winner The array to store the winner variables
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[idx].has_shower_ended()) return;
  event& ev = events[idx];
  // ---------------------------------------------

  // default values
  double win_tt = shower->t_c;  // lowest possible value is cutoff
  int win_sf = 0;               // 0 = no splitting
  int win_ij = 0;
  int win_k = 0;
  double win_sijk = 0.;
  double win_zm = 0.;
  double win_zp = 0.;

  for (int ij = 0; ij < ev.get_size(); ij++) {
    for (int k = 0; k < ev.get_size(); k++) {
      // sanity check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // Skip non-partons
      if (!ev.get_particle(ij).is_parton() || !ev.get_particle(k).is_parton()) {
        continue;
      }

      // need to check if ij and k are colour connected
      if (!ev.get_particle(ij).is_color_connected(ev.get_particle(k))) {
        continue;
      }

      // get the invariant mass squared of the dipole
      double sijk =
          (ev.get_particle(ij).get_mom() + ev.get_particle(k).get_mom()).m2();

      // get the splitting functions for the current partons
      int sf_codes[11];
      shower->generate_possible_splittings(
          ev.get_particle(ij).get_pid(), ev.get_particle(k).get_pid(),
          ev.get_particle(ij).is_initial(), ev.get_particle(k).is_initial(),
          sf_codes);

      // codes instead of object oriented approach!
      for (int sf : sf_codes) {
        // When a null code is encountered, we have reached the end of the
        // possible splittings, we can break out of the loop
        if (sf == -1) {
          break;
        }

        // check if either particle has eta < eta_min (usually 1e-5)
        if (shower->is_ii(sf) && (ev.get_particle(ij).get_eta() < 1e-5 ||
                                  ev.get_particle(k).get_eta() < 1e-5)) {
          continue;
        }

        // phase space limits
        double zm, zp;
        double eta = shower->is_fi(sf) ? ev.get_particle(k).get_eta()
                                       : ev.get_particle(ij).get_eta();
        shower->get_boundaries(zm, zp, sijk, eta, sf);
        if (zm < 0. || zp > 1. || zm > zp) {
          continue;
        }

        // calculate the integrated overestimate
        double pdf_max = shower->get_pdf_max(sf, ev.get_particle(ij).get_eta());
        double c = shower->as_max / (2. * M_PI) *
                   shower->sf_integral(zm, zp, sf) * shower->j0_max * pdf_max;

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

          tt = (ev.get_shower_t() - mq2) * pow(ev.gen_random(), 1. / c) + mq2;

          // Check if tt <= mq2: Need some numerical error protection
          if (tt - mq2 < 1e-9) {
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

__global__ void check_cutoff(event* events, shower* shower, int* d_completed,
                             int n) {
  /**
   * @brief Check if the shower has ended
   *
   * @param events The events to run the shower on
   * @param d_completed The number of completed events
   * @param cutoff The cutoff scale
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[idx].has_shower_ended()) return;
  event& ev = events[idx];
  // ---------------------------------------------

  /**
   * end shower if t < cutoff
   *
   * ev.get_shower_t() <= cutoff is equally valid
   * i just prefer this way because this way is
   * how we usually write it in literature
   */
  if (!(ev.get_shower_t() > shower->t_c)) {
    ev.shower_has_ended(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events

    return;
  }
}

// -----------------------------------------------------------------------------

/**
 * PDF Ratio Calculation
 * ---------------------
 *
 * This is done in two steps:
 * - PDFs are evaluated for ij and i, see pdf.cuh
 * - Ratio is calculated in the veto_alg kernel
 */

// -----------------------------------------------------------------------------

__global__ void veto_alg(shower* shower, alpha_s* as, event* events, int n,
                         double* xf_a, double* xf_b, bool* accept_emission,
                         double* winner) {
  /**
   * @brief The veto algorithm for the shower
   *
   * @param shower The shower object
   * @param as The alpha_s object
   * @param events The events to run the shower on
   * @param n The number of events
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
  if (events[idx].has_shower_ended()) return;
  event& ev = events[idx];
  // ---------------------------------------------

  // set to false, only set to true if accpeted
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

  // Check Phase Space is Valid
  if (!(shower->check_phase_space(z, y, sf))) {
    return;
  }

  // Additional Check for FI and IF/II, in case quark > proton
  vec4 pijt = ev.get_particle(ij).get_mom();
  vec4 pkt = ev.get_particle(k).get_mom();
  if ((shower->is_isr(sf)) && (pijt[0] / z > shower->e_proton)) {
    return;
  }
  if ((shower->is_fi(sf)) && (pkt[0] / y > shower->e_proton)) {
    return;
  }

  // Get PDF Ratio and PDF Max for FI and IF/II
  double pdf_ratio(1.), pdf_max(1.);
  if (!shower->is_ff(sf)) {
    // Check Momentum Fraction
    if (!(shower->check_mom_frac(
            sf, ev.get_particle(ij).get_pid(), ev.get_particle(k).get_pid(),
            ev.get_particle(ij).get_eta(), ev.get_particle(k).get_eta(), z))) {
      return;
    }

    // Ensure PDFs are valid
    if (isnan(xf_a[idx]) || isnan(xf_b[idx]) || isinf(xf_a[idx]) ||
        isinf(xf_b[idx]) || xf_a[idx] <= 0. || xf_b[idx] <= 0.) {
      return;
    }

    // Calculate the ratio and handle division by zero
    pdf_ratio = xf_a[idx] / xf_b[idx];

    // Verify the pdfratio
    if (isnan(pdf_ratio) || isinf(pdf_ratio) || (pdf_ratio <= 0.)) {
      return;
    }

    // If PDF Ratio is too large, veto
    if (pdf_ratio > 1e6) {
      return;
    }

    // Get PDF Max
    pdf_max = shower->get_pdf_max(sf, ev.get_particle(ij).get_eta());

    /**
     * Why the factor of z?
     * --------------------
     * From LHAPDF, we get xf(x, q2). This means our ratio will be equal to
     * (x/z)f(x/z, q2) / xf(x, q2) = 1/z * f(x/z, q2) / f(x, q2)
     *
     * So, we can either adjust the jacobian by multiplying by z or adjust the
     * pdf_ratio by dividing by z. We choose the latter.
     */
    pdf_ratio *= z;

    // Mutliply by (t - m2) / t for ISR to account for quark masses
    int fl = abs(ev.get_particle(ij).get_pid());
    pdf_ratio *= (t - (fl == 5 ? mb * mb : (fl == 4 ? mc * mc : 0.))) / t;
  }

  // Jacobian
  double jacobian = shower->get_jacobian(z, y, sf) * pdf_ratio;
  double jmaxtot = shower->j0_max * pdf_max;

  // Splitting Function Value and Estimate
  double value = shower->sf_value(z, y, sf);
  double estimate = shower->sf_estimate(z0, sf);

  // veto algorithm
  double f = (*as)(t)*value * jacobian;
  double g = shower->as_max * estimate * jmaxtot;

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
__global__ void do_splitting(shower* shower, event* events, int n,
                             bool* accept_emission, double* winner) {
  /**
   * @brief Do the splitting for the shower
   *
   * @param shower The shower object
   * @param events The events to run the shower on
   * @param n The number of events
   * @param accept_emission The array to store the acceptance of the emission
   * @param winner The array to store the winner emission data
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[idx].has_shower_ended()) return;
  event& ev = events[idx];
  // ---------------------------------------------

  // Do not run if the shower has ended
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

  // get flavs from kernel number
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

  // add emitted particle
  particle em = particle(flavs[2], moms[1], colj[0], colj[1]);
  ev.set_particle(ev.get_size(), em);

  // increment emissions (important) !!!!!
  ev.increment_emissions();

  // II Only - Lorentz Boost the new final state
  if (shower->is_ii(sf)) {
    shower->ii_boost_after_emission(ev, moms);
  }

  return;
}

// -----------------------------------------------------------------------------

__global__ void check_too_many_particles(event* events, int n_emissions_max,
                                         int* d_too_many_particles,
                                         int* d_completed, int n) {
  /**
   * @brief Check if the event has too many particles
   *
   * @param events The events to run the shower on
   * @param d_too_many_particles The number of events with too many particles
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Shower Preamble
  if (events[idx].has_shower_ended()) return;
  event& ev = events[idx];
  // ---------------------------------------------

  // limit to max particles
  if (ev.get_size() == min(max_particles, ev.get_hard() + n_emissions_max)) {
    ev.shower_has_ended(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events
    atomicAdd(d_too_many_particles, 1);  // to print later
    return;
  }
}

// -----------------------------------------------------------------------------

struct is_not_end_shower {
  /**
   * @brief Function object to check if the shower has ended
   *
   * @param ev The event to check
   * @return true if the shower has not ended
   */
  __device__ bool operator()(const event& ev) const {
    return !ev.has_shower_ended();
  }
};

// -----------------------------------------------------------------------------

void run_shower(thrust::device_vector<event>& dv_events, double root_s,
                bool nlo_matching, bool do_partitioning, double t_c,
                double asmz, bool fixed_as, int n_emissions_max, int blocks,
                int threads) {
  /**
   * @brief Run the shower on the events
   *
   * @param dv_events The events to run the shower on
   * @param root_s The root s energy
   * @param nlo_matching Whether to do NLO matching
   * @param do_partitioning Whether to partition the events
   */

  // number of events - can get from d_events.size()
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // set up the device alpha_s calculator
  alpha_s* d_as;
  cudaMalloc(&d_as, sizeof(alpha_s));
  as_setup_kernel<<<1, 1>>>(d_as, mz, asmz, (fixed_as ? 0 : 2));
  sync_gpu_and_check("as_setup_kernel");

  // Calculate as_max = as(t_c)
  double* d_as_max;
  cudaMalloc(&d_as_max, sizeof(double));
  as_value<<<1, 1>>>(d_as, d_as_max, t_c);
  sync_gpu_and_check("as_value");
  double as_max;
  cudaMemcpy(&as_max, d_as_max, sizeof(double), cudaMemcpyDeviceToHost);

  // set up the shower
  shower* d_shower;
  cudaMalloc(&d_shower, sizeof(shower));
  shower_setup_kernel<<<1, 1>>>(d_shower, root_s, t_c, as_max);
  sync_gpu_and_check("shower_setup_kernel");

  // set up the pdf evaluator
  pdf_wrapper pdf;

  /**
   * Shower Variables - useful to store as collective
   *
   * t, c and end_shower stored in event, becuase they
   * are unique to each event, and not throwaway values
   * like these.
   *
   * Winner variables: (sf, ij, k, sijk, z, y, phi)
   * Stored in ONE array, so we make it 6 x n_events
   * Stored all as doubles, so static_cast<int> for sf, ij, k
   */
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

  // ---------------------------------------------------
  // Analysis Variables

  // allocate device memory for completed events counter
  int* d_completed;
  cudaMalloc(&d_completed, sizeof(int));
  cudaMemset(d_completed, 0, sizeof(int));

  // allocate device memory to counts events that surpass max particles
  int* d_too_many_particles;
  cudaMalloc(&d_too_many_particles, sizeof(int));
  cudaMemset(d_too_many_particles, 0, sizeof(int));

  // ---------------------------------------------------------------------------
  // prepare the shower

  debug_msg("running @prep_shower");
  prep_shower<<<blocks, threads>>>(d_events, nlo_matching, n_events);
  sync_gpu_and_check("prep_shower");

  // ---------------------------------------------------------------------------
  // run the shower

  // number of completed events and cycles
  int completed = 0;
  int cycle = 0;

  // (Varying) kernel size and partition factor
  int n = n_events;
  int p = 1;

  while (completed < n_events) {
    // run all the kernels here...

    // -------------------------------------------------------------------------
    // check if there are too many particles (do first in case of H event)

    debug_msg("running @check_too_many_particles");
    check_too_many_particles<<<blocks, threads>>>(
        d_events, n_emissions_max, d_too_many_particles, d_completed, n);
    sync_gpu_and_check("check_too_many_particles");

    // -------------------------------------------------------------------------
    // select the winner kernel

    debug_msg("running @select_winner_split_func");
    select_winner_split_func<<<blocks, threads>>>(d_shower, d_events, n,
                                                  d_winner);
    sync_gpu_and_check("select_winner_split_func");

    // -------------------------------------------------------------------------
    // check cutoff

    debug_msg("running @check_cutoff");
    check_cutoff<<<blocks, threads>>>(d_events, d_shower, d_completed, n);
    sync_gpu_and_check("check_cutoff");

    // -------------------------------------------------------------------------
    // calculate pdf of ij and i using LHAPDF

    debug_msg("running @setup_pdfratio");
    setup_pdfratio<<<blocks, threads>>>(d_shower, d_events, n, d_flavours_a,
                                        d_flavours_b, d_x_a, d_x_b, d_q2,
                                        d_winner);
    sync_gpu_and_check("setup_pdfratio");

    pdf.evaluate(d_flavours_a, d_x_a, d_q2, d_xf_a, n, blocks, threads);
    pdf.evaluate(d_flavours_b, d_x_b, d_q2, d_xf_b, n, blocks, threads);

    // -------------------------------------------------------------------------
    // veto algorithm

    debug_msg("running @veto_alg");
    veto_alg<<<blocks, threads>>>(d_shower, d_as, d_events, n, d_xf_a, d_xf_b,
                                  d_accept_emission, d_winner);
    sync_gpu_and_check("veto_alg");

    // -------------------------------------------------------------------------
    // splitting algorithm

    debug_msg("running @do_splitting");
    do_splitting<<<blocks, threads>>>(d_shower, d_events, n, d_accept_emission,
                                      d_winner);
    sync_gpu_and_check("do_splitting");

    // -------------------------------------------------------------------------
    // import the number of completed events

    cudaMemcpy(&completed, d_completed, sizeof(int), cudaMemcpyDeviceToHost);
    cycle++;

    // -------------------------------------------------------------------------

    // Partition the events based on completion at 50%, 75%, 87.5%, etc.
    if (do_partitioning &&
        completed >= static_cast<int>(n_events * (1 - 1 / pow(2, p)))) {
      if (static_cast<int>(n_events * (1 / pow(2, p))) > 10) {
        std::cerr << std::endl;
        std::cerr << "partition at " << completed << "/" << n_events
                  << std::endl;

        // Partition only the first n elements
        // Prevents partitioning of completed events
        // using n ensures all incomplete events are
        // accounted for
        thrust::partition(dv_events.begin(), dv_events.begin() + n,
                          is_not_end_shower());

        // Update n to reflect the number of incomplete events
        n = n_events - completed;

        // Increment the partitioning step
        p++;
      }
    }

    // -------------------------------------------------------------------------
    // show progress

    std::cerr << "\rCompleted Events: " << completed << "/" << n_events
              << std::flush;

    // -------------------------------------------------------------------------
  }
  std::cout << std::endl;

  // print the number of events that surpassed the max particles
  int too_many_particles;
  cudaMemcpy(&too_many_particles, d_too_many_particles, sizeof(int),
             cudaMemcpyDeviceToHost);
  if (too_many_particles > 0) {
    if (max_particles < n_emissions_max) {
      std::cerr << "Warning: " << too_many_particles
                << " events surpassed the maximum number of particles"
                << std::endl;
      std::cerr << "Consider increasing max_particles, default: "
                << max_particles << std::endl;
    }
  }

  // ---------------------------------------------------------------------------

  // free the memory
  cudaFree(d_shower);
  cudaFree(d_as);
  cudaFree(d_completed);
  cudaFree(d_too_many_particles);
}