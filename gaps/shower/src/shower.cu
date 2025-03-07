#include "shower.cuh"

// Needed for Checking Time Taken for Cycles
// #include <chrono>
// #include <fstream>

// need to be here to avoid multiple definitions
#include "colours.cuh"
#include "kinematics.cuh"
#include "splittings.cuh"

// -----------------------------------------------------------------------------
// preparing the shower

__global__ void prep_shower(event *events, int n) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if (idx >= n) {
    return;
  }

  event &ev = events[idx];

  // set the starting shower scale
  double t_max = (ev.get_parton(0).get_mom() + ev.get_parton(1).get_mom()).m2();
  ev.set_shower_t(t_max);

  // set the initial number of emissions
  ev.set_emissions(0);

  // set the colour counter to 1 (q and qbar)
  ev.set_shower_c(1);

  // set the initial end shower flag - no need, default value is false
  // ev.set_end_shower(false);
}

// -----------------------------------------------------------------------------

/**
 * selecting the winner splitting function
 * ---------------------------------------
 *
 * when you profile the code, you will notice that this is the process that
 * takes up half of the shower time. this method below is a first attempt at
 * parallelizing the process.
 */

__global__ void select_winner_split_func(event *events, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event &ev = events[idx];

  // do not run if the shower has ended
  if (ev.get_end_shower()) {
    return;
  }

  // default values
  double win_tt = t_c;  // lowest possible value is cutoff scale (in base.cuh)
  int win_sf = 0;       // 0 = no splitting
  int win_ij = 0;
  int win_k = 0;
  double win_zp = 0.;
  double win_m2 = 0.;

  // we start at 2 because elements 0 and 1 are electrons - to change with isr
  for (int ij = 2; ij < ev.get_size(); ij++) {
    for (int k = 2; k < ev.get_size(); k++) {
      // sanity check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // need to check if ij and k are colour connected
      if (!ev.get_parton(ij).is_color_connected(ev.get_parton(k))) {
        continue;
      }

      // params identical to all splitting functions
      double m2 =
          (ev.get_parton(ij).get_mom() + ev.get_parton(k).get_mom()).m2();
      if (m2 < 4. * t_c) {
        continue;
      }

      // phase space limits
      double zp = 0.5 * (1. + sqrt(1. - 4. * t_c / m2));

      // generate the splitting function code:
      int sf_codes[6] = {-1, -1, -1, -1, -1, -1};
      if (ev.get_parton(ij).get_pid() != 21) {
        // quark
        if (ev.get_parton(ij).get_pid() > 0) {
          sf_codes[0] = ev.get_parton(ij).get_pid();
        }
        // antiquark
        else {
          // -1 -> 11, -2 -> 12, -3 -> 13, -4 -> 14, -5 -> 15
          sf_codes[0] = -ev.get_parton(ij).get_pid() + 10;
        }
      }
      // gluon
      else {
        sf_codes[0] = 200;
        sf_codes[1] = 301;
        sf_codes[2] = 302;
        sf_codes[3] = 303;
        sf_codes[4] = 304;
        sf_codes[5] = 305;
      }

      // codes instead of object oriented approach!
      for (int sf : sf_codes) {
        // check if the splitting function is valid for the current partons
        if (sf == -1) {
          break;
        }

        // calculate the evolution variable
        double g = asmax / (2. * M_PI) * sf_integral(1 - zp, zp, sf);
        double tt = ev.get_shower_t() * pow(ev.gen_random(), 1. / g);

        // check if tt is greater than the current winner
        if (tt > win_tt) {
          win_tt = tt;
          win_sf = sf;
          win_ij = ij;
          win_k = k;
          win_zp = zp;
          win_m2 = m2;
        }
      }
    }
  }

  // store the results
  ev.set_shower_t(win_tt);
  ev.set_win_sf(win_sf);
  ev.set_win_dipole(0, win_ij);
  ev.set_win_dipole(1, win_k);
  ev.set_win_param(0, win_zp);
  ev.set_win_param(1, win_m2);
}

// -----------------------------------------------------------------------------

__global__ void check_cutoff(event *events, int *d_completed, double cutoff,
                             int *d_reached_max_partons, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event &ev = events[idx];

  // do not run if the shower has ended
  if (ev.get_end_shower()) {
    return;
  }

  // limit to max_partons
  if (ev.get_size() == max_partons) {
    ev.set_end_shower(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events
    atomicAdd(d_reached_max_partons, 1);  // increment the number of events that
                                           // reached the maximum number of partons
    return;
  }
  

  /**
   * end shower if t < cutoff
   *
   * ev.get_shower_t() <= cutoff is equally valid
   * i just prefer this way because this way is
   * how we usually write it in literature
   */
  if (!(ev.get_shower_t() > cutoff)) {
    ev.set_end_shower(true);
    atomicAdd(d_completed, 1);  // increment the number of completed events

    return;
  }
}

// -----------------------------------------------------------------------------

__global__ void veto_alg(event *events, double *asval, bool *accept_emission,
                         int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event &ev = events[idx];

  // do not run if the shower has ended
  if (ev.get_end_shower()) {
    return;
  }

  // set to false, only set to true if accpeted
  accept_emission[idx] = false;

  // get the splitting function
  int sf = ev.get_win_sf();

  // generate z
  double zp = ev.get_win_param(0);
  double z = sf_generate_z(1 - zp, zp, ev.gen_random(), sf);

  double y = ev.get_shower_t() / ev.get_win_param(1) / z / (1. - z);

  double f = 0.;
  double g = 0.;
  double value = 0.;
  double estimate = 0.;

  // cs kernel: y can't be 1
  if (y < 1.) {
    value = sf_value(z, y, sf);
    estimate = sf_estimate(z, sf);

    f = (1. - y) * asval[idx] * value;
    g = asmax * estimate;

    if (ev.gen_random() < f / g) {
      accept_emission[idx] = true;
      ev.set_shower_z(z);
      ev.set_shower_y(y);
    }
  }
}

// -----------------------------------------------------------------------------

// do splitting
__global__ void do_splitting(event *events, bool *accept_emission, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event &ev = events[idx];

  // do not run if the shower has ended
  if (ev.get_end_shower()) {
    return;
  }

  if (!accept_emission[idx]) {
    return;
  }

  double phi = 2. * M_PI * ev.gen_random();

  int win_ij = ev.get_win_dipole(0);
  int win_k = ev.get_win_dipole(1);

  // make kinematics
  vec4 moms[3] = {vec4(), vec4(), vec4()};

  make_kinematics(moms, ev.get_shower_z(), ev.get_shower_y(), phi,
                  ev.get_parton(win_ij).get_mom(),
                  ev.get_parton(win_k).get_mom());

  // adjust colors

  // get flavs from kernel number
  int sf = ev.get_win_sf();
  int flavs[3];
  sf_to_flavs(sf, flavs);

  int colij[2] = {ev.get_parton(win_ij).get_col(),
                  ev.get_parton(win_ij).get_anti_col()};

  int colk[2] = {ev.get_parton(win_k).get_col(),
                 ev.get_parton(win_k).get_anti_col()};

  int coli[2] = {0, 0};
  int colj[2] = {0, 0};

  make_colours(ev, coli, colj, flavs, colij, colk, ev.gen_random());

  // modify splitter
  ev.set_parton_pid(win_ij, flavs[1]);
  ev.set_parton_mom(win_ij, moms[0]);
  ev.set_parton_col(win_ij, coli[0]);
  ev.set_parton_anti_col(win_ij, coli[1]);

  // modify recoiled spectator
  ev.set_parton_mom(win_k, moms[2]);

  // add emitted parton
  parton em = parton(flavs[2], moms[1], colj[0], colj[1]);
  ev.set_parton(ev.get_size(), em);

  // increment emissions (important)
  ev.increment_emissions();

  return;
}

// -----------------------------------------------------------------------------

struct is_not_end_shower {
  /**
   * @brief Function object to check if the shower has ended
   *
   * @param ev The event to check
   * @return true if the shower has not ended
   */
  __device__ bool operator()(const event &ev) const {
    return !ev.get_end_shower();
  }
};

// -----------------------------------------------------------------------------

void run_shower(thrust::device_vector<event> &dv_events, bool do_partition) {
  // number of events - can get from d_events.size()
  event *d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  int threads = 128;
  int blocks = (n_events + threads - 1) / threads;

  // set up the device alpha_s
  alpha_s *d_as;
  cudaMalloc(&d_as, sizeof(alpha_s));
  as_setup_kernel<<<1, 1>>>(d_as, mz, asmz);
  sync_gpu_and_check("as_setup_kernel");

  // allocate device memory for completed events counter
  int *d_completed;
  cudaMalloc(&d_completed, sizeof(int));
  cudaMemset(d_completed, 0, sizeof(int));

  // allocate device memory for the number of events that reached the maximum
  // number of partons
  int *d_reached_max_partons;
  cudaMalloc(&d_reached_max_partons, sizeof(int));
  cudaMemset(d_reached_max_partons, 0, sizeof(int));

  // as(t) and veto
  double *d_asval;
  cudaMalloc(&d_asval, n_events * sizeof(double));
  bool *d_accept_emission;
  cudaMalloc(&d_accept_emission, n_events * sizeof(bool));

  // store the number of finished events per cycle
  std::vector<int> completed_per_cycle;


  // ---------------------------------------------------------------------------
  // prepare the shower

  debug_msg("running @prep_shower");
  prep_shower<<<blocks, threads>>>(d_events, n_events);
  sync_gpu_and_check("prep_shower");

  // ---------------------------------------------------------------------------
  // run the shower

  int completed = 0;
  int cycle = 0;
  // auto start = std::chrono::high_resolution_clock::now();
  // auto end = std::chrono::high_resolution_clock::now();
  // double diff = 0.;

  // (Varying) kernel size and partition factor
  int n = n_events;
  int p = 1;

  while (completed < n_events) {
    // run all the kernels here...

    // start the clock
    // start = std::chrono::high_resolution_clock::now();

    // -------------------------------------------------------------------------
    // select the winner kernel

    debug_msg("running @select_winner_split_func");
    select_winner_split_func<<<blocks, threads>>>(d_events, n);
    sync_gpu_and_check("select_winner_split_func");

    // -------------------------------------------------------------------------
    // check cutoff

    debug_msg("running @check_cutoff");
    check_cutoff<<<blocks, threads>>>(d_events, d_completed, t_c, 
                                            d_reached_max_partons, n);
    sync_gpu_and_check("check_cutoff");

    // -------------------------------------------------------------------------
    // calculate alpha_s for veto algorithm

    debug_msg("running @as_shower_kernel");
    as_shower_kernel<<<blocks, threads>>>(d_as, d_events, d_asval, n);
    sync_gpu_and_check("as_shower_kernel");

    // -------------------------------------------------------------------------
    // veto algorithm

    debug_msg("running @veto_alg");
    veto_alg<<<blocks, threads>>>(d_events, d_asval, d_accept_emission,
                                       n);
    sync_gpu_and_check("veto_alg");

    // -------------------------------------------------------------------------
    // splitting algorithm

    debug_msg("running @do_splitting");
    do_splitting<<<blocks, threads>>>(d_events, d_accept_emission, n);
    sync_gpu_and_check("do_splitting");

    // -------------------------------------------------------------------------
    // import the number of completed events

    cudaMemcpy(&completed, d_completed, sizeof(int), cudaMemcpyDeviceToHost);
    cycle++;

    // end the clock
    // end = std::chrono::high_resolution_clock::now();
    // diff = std::chrono::duration<double>(end - start).count();
    // time_per_cycle.push_back(diff);

    // until paper is published, we will use this
    completed_per_cycle.push_back(completed);
    std::cerr << "\rCompleted Events: " << completed << "/" << n_events
    << std::flush;

    // -------------------------------------------------------------------------
   
    // Partition the events based on completion at 50%, 75%, 87.5%, etc.
    if (do_partition &&
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
  }
  std::cerr << std::endl;
  // ---------------------------------------------------------------------------
  // write completed_per_cycle to file
  std::ofstream file("gaps-cycles.dat");
  for (auto &i : completed_per_cycle) {
    file << i << std::endl;
  }

  int h_reached_max_partons = 0;
  cudaMemcpy(&h_reached_max_partons, d_reached_max_partons, sizeof(int),
             cudaMemcpyDeviceToHost);
  std::cout << "number of events that reached the maximum number of partons: "
            << h_reached_max_partons << std::endl;

  // for (int i = 0; i < cycle; i++) {
  //   file << i << ", " << completed_per_cycle[i] << ", " << time_per_cycle[i]
  //        << std::endl;
  // }

  // ---------------------------------------------------------------------------
  // clean up device memory
  cudaFree(d_asval);
  cudaFree(d_accept_emission);
  cudaFree(d_completed);
}
