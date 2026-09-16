#include "hadronic_decays.cuh"

// -----------------------------------------------------------------------------
// Stability Check

__device__ bool hadronic_decays::is_stable(int pid) const {
  /**
   * @brief Is this a final state particle?
   *
   * @param pid the particle id
   */

  // Unknown ids are final; a state with no width (K0, K0bar) always decays
  int i = hadron_index(pid);
  if (i < 0) return true;
  double ctau = hadron_database[i].ctau;
  if (ctau == ctau_no_width) return false;
  return ctau >= ctau_max;
}

// -----------------------------------------------------------------------------
// Single Decay

__device__ bool hadronic_decays::decay_one(event& ev, int i) const {
  /**
   * @brief Decay the state at index i into its children.
   *
   * @param ev the event to decay in
   * @param i the index of the state to decay
   * @return false where the state was left as it was
   */

  // Parent
  particle parent = ev.get_particle(i);

  // A state the table does not hold, or one with no mode of its own, is left
  // standing
  int row = hadron_index(parent.get_pid());
  if (row < 0 || hadron_database[row].n_modes == 0) return false;

  // Its modes, which the table keeps together
  const decay_row* modes = decay_database + hadron_database[row].mode_begin;
  int n_modes = hadron_database[row].n_modes;

  // Throw one mode with probability br / sum(br), as choose_with_weights does
  // but walking the rows in place, so the table need not be normalised
  double total = 0.;
  for (int m = 0; m < n_modes; m++) total += modes[m].br;

  double rho_mode = ev.gen_random();
  double target = rho_mode * total;

  int chosen = n_modes - 1;
  double running = 0.;
  for (int m = 0; m < n_modes; m++) {
    running += modes[m].br;
    if (running > target) {
      chosen = m;
      break;
    }
  }

  const decay_row& mode = modes[chosen];
  int n_children = mode.n_children;

  // The K0 and the K0bar mix into one child, which takes the whole momentum
  // under a new id
  if (n_children == 1) {
    ev.set_particle(i, particle(mode.children[0], parent.get_mom(), 0, 0));
    return true;
  }

  // The masses the children carry, every child being a state of the table
  double masses[max_decay_products];
  for (int c = 0; c < n_children; c++) {
    masses[c] = hadron_mass(mode.children[c]);
  }

  // Run the decay kinematics, one chain of two-body throws
  vec4 momenta[max_decay_products];
  if (!one_to_n_decay(parent.get_mom(), masses, n_children, momenta, ev)) {
    printf("hadronic_decays: %d of mass %f GeV left undecayed\n",
           parent.get_pid(), parent.get_mom().m());
    return false;
  }

  // Replace the parent, and append the rest
  ev.set_particle(i, particle(mode.children[0], momenta[0], 0, 0));

  for (int c = 1; c < n_children; c++) {
    ev.add_particle(particle(mode.children[c], momenta[c], 0, 0));
  }

  return true;
}

// -----------------------------------------------------------------------------
// Decay Loop

__device__ void hadronic_decays::run(event& ev) const {
  /**
   * @brief Decay every unstable particle until the event is final.
   *
   * @param ev the event to decay
   */

  for (int generation = 0; generation < max_generations; generation++) {
    bool all_stable = true;

    // children appended below wait for the next pass
    int n_particles = ev.get_size();

    // Check if stable, or decay
    for (int i = 2; i < n_particles; i++) {
      if (is_stable(ev.get_particle(i).get_pid())) continue;

      if (decay_one(ev, i)) all_stable = false;
    }

    if (all_stable) return;
  }

  printf("hadronic_decays: event still holds unstable particles after %d "
         "generations\n",
         max_generations);
}

// -----------------------------------------------------------------------------
// Kernel

__global__ void hadronic_decays_kernel(event* events, hadronic_decays* dec,
                                       int n) {
  /**
   * @brief Decay the hadrons of every event.
   *
   * @param events the events to decay
   * @param dec the decay object
   * @param n the number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  // Check for overflow
  if (ev.get_overflowed()) return;
  // ---------------------------------------------

  dec->run(ev);
}

// -----------------------------------------------------------------------------
// Wrapper

void run_decays(thrust::device_vector<event>& dv_events, const params& p,
                int blocks) {
  /**
   * @brief Run the hadron decays on the events
   *
   * @param dv_events The device vector of hadronised events
   * @param p         Run parameters (threads, ctau_max, ...)
   * @param blocks    Number of CUDA blocks
   */

  // Get raw pointer and event count from the thrust vector
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // Set up the device decay object
  hadronic_decays *h_dec, *d_dec;
  h_dec = new hadronic_decays(p.ctau_max);
  cudaMalloc(&d_dec, sizeof(hadronic_decays));
  cudaMemcpy(d_dec, h_dec, sizeof(hadronic_decays), cudaMemcpyHostToDevice);

  // Decay everything the cluster decay left unstable
  debug_msg("running @hadronic_decays");
  hadronic_decays_kernel<<<blocks, p.threads>>>(d_events, d_dec, n_events);
  sync_gpu_and_check("hadronic_decays");

  // Clean up
  delete h_dec;
  cudaFree(d_dec);
}
