#include "shower.cuh"

// -----------------------------------------------------------------------------
// Random Number Generator

// No need during matrix as initialised once and used once only
// But for shower used 80 to 100 times
__global__ void initCurandStates(curandState *states, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= N) {
    return;
  }
  // Every events[idx] has a seed idx
  //curand_init(idx, 0, 0, &states[idx]);

  // Every events[idx] has a seed idx and clok64() is used to get a seed
  curand_init(clock64(), idx, 0, &states[idx]);
}

// -----------------------------------------------------------------------------
// Preparing the Shower

__global__ void prepShower(Event *events, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  if (idx >= N) {
    return;
  }

  Event &ev = events[idx];

  // Set the starting shower scale
  double t_max = (ev.GetParton(0).GetMom() + ev.GetParton(1).GetMom()).M2();
  ev.SetShowerT(t_max);

  // Set the initial number of emissions
  ev.SetEmissions(0);

  // Set the Colour Counter to 1 (q and qbar)
  ev.SetShowerC(1);

  // Set the initial end shower flag - No need, default value is false
  // ev.SetEndShower(false);
}

// -----------------------------------------------------------------------------

/**
 * Selecting the Winner Splitting Function
 * ---------------------------------------
 *
 * When you profile the code, you will notice that this is the process that
 * takes up half of the shower time. This method below is a first attempt at
 * parallelizing the process.
 */

__global__ void selectWinnerSplitFunc(Event *events, curandState *states,
                                      int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  curandState state = states[idx];

  Event &ev = events[idx];

  // Do not run if the shower has ended
  if (ev.GetEndShower()) {
    return;
  }

  // Default Values
  double win_tt = tC;  // Lowest possible value is Cutoff Scale (in base.cuh)
  int win_sf = 16;     // 16 = No Splitting (0 -> 15 are Splitting Functions)
  int win_ij = 0;
  int win_k = 0;
  double win_zp = 0.0;
  double win_m2 = 0.0;

  for (int ij = 2; ij < ev.GetSize(); ij++) {
    for (int k = 2; k < ev.GetSize(); k++) {
      // Sanity Check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // Need to check if ij and k are valid partons
      if (ij >= ev.GetSize() || k >= ev.GetSize()) {
        continue;
      }

      // Need to check if ij and k are colour connected
      if (!ev.GetParton(ij).IsColorConnected(ev.GetParton(k))) {
        continue;
      }

      // Params Identical to all splitting functions
      double m2 = (ev.GetParton(ij).GetMom() + ev.GetParton(k).GetMom()).M2();
      if (m2 < 4.0 * tC) {
        continue;
      }

      double zp = 0.5 * (1.0 + sqrt(1.0 - 4.0 * tC / m2));

      /**
       * SIMPLIFIED VERSION FROM S. HOCHE'S CODE
       * ---------------------------------------
       *
       * CUDA doesn't directly support Polymorphism, so we can't use the
       * SplittingFunction class and its subclasses here. We will have to
       * implement the splitting functions directly here...
       *
       * There might be possible workarounds for this:
       * - Use a template function
       * - Use nvstd::function
       *
       * For now (Feb 2024), let's keep it simple and implement the splitting
       * functions directly here. It is good that we are doing massless quarks.
       *
       * HOWEVER this change is present in both the C++ and CUDA showers, so
       * the comparison is still valid.
       *
       * If ij = quark/antiq, only one of the first 10 kernels can be used
       * If ij = gluon, we need to get one of six kernels
       */

      int sf = 16;  // Default or Null Splitting Function
      double tt = 0.0;

      // Quark/Anti-Quark
      if (ev.GetParton(ij).GetPid() != 21) {
        // Flavour:  -5 -4 -3 -2 -1  1  2  3  4  5
        // Location:  0  1  2  3  4  5  6  7  8  9

        // Get the right splitting function
        sf = ev.GetParton(ij).GetPid() < 0 ? ev.GetParton(ij).GetPid() + 5
                                           : ev.GetParton(ij).GetPid() + 4;

        // Get t_temp
        double integral = kCF * 2.0 * log((1.0 - (1.0 - zp)) / (1.0 - zp));
        double g = asmax / (2.0 * M_PI) * integral;
        tt = ev.GetShowerT() * pow(curand_uniform(&state), 1.0 / g);
        states[idx] = state;  // So that the next number is not the same!
      } else if (ev.GetParton(ij).GetPid() == 21) {
        /**
         * We will need to give a chance to all six possible splittings
         *
         * 10: g -> gg
         * 11 to 15: g -> qqbar
         */

        double gg_integral = kCA * log((1.0 - (1.0 - zp)) / (1.0 - zp));
        double qq_integral = kTR / 2.0 * (zp - (1.0 - zp));

        double gg_g = asmax / (2.0 * M_PI) * gg_integral;
        double qq_g = asmax / (2.0 * M_PI) * qq_integral;

        // Get tt for g -> gg
        double gg_tt =
            ev.GetShowerT() * pow(curand_uniform(&state), 1.0 / gg_g);
        states[idx] = state;

        // Get tt for g -> qqbar
        // Select the quark flavour with highest random number
        double qq_tt = 0.0;
        int qnum = 0;
        for (int j = 1; j <= 5; j++) {
          // Get tt for g -> qqbar for the selected quark flavour
          double temp_tt =
              ev.GetShowerT() * pow(curand_uniform(&state), 1.0 / qq_g);
          states[idx] = state;

          if (temp_tt > qq_tt) {
            qq_tt = temp_tt;
            qnum = j;
          }
        }

        // Compare tt for g -> gg and g -> qqbar
        if (gg_tt > qq_tt) {
          sf = 10;
          tt = gg_tt;
        } else {
          sf = 10 + qnum;
          tt = qq_tt;
        }
      }

      // Check if tt is greater than the current winner
      if (sf < 16 && tt > win_tt) {
        win_tt = tt;
        win_sf = sf;
        win_ij = ij;
        win_k = k;
        win_zp = zp;
        win_m2 = m2;
      }
    }
  }

  // Store the results
  ev.SetShowerT(win_tt);
  ev.SetWinSF(win_sf);
  ev.SetWinDipole(0, win_ij);
  ev.SetWinDipole(1, win_k);
  ev.SetWinParam(0, win_zp);
  ev.SetWinParam(1, win_m2);
}

// -----------------------------------------------------------------------------

__global__ void checkCutoff(Event *events, int *d_completed, double cutoff,
                            int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event &ev = events[idx];

  // Do not run if the shower has ended
  if (ev.GetEndShower()) {
    return;
  }

  /**
   * End shower if t < cutoff
   *
   * ev.GetShowerT() <= cutoff is equally valid
   * I just prefer this way because this way is
   * how we usually write it in Literature
   */
  if (!(ev.GetShowerT() > cutoff)) {
    ev.SetEndShower(true);
    atomicAdd(d_completed, 1);  // Increment the number of completed events
  }
}

// -----------------------------------------------------------------------------

__global__ void vetoAlg(Event *events, curandState *states, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event &ev = events[idx];

  // Do not run if the shower has ended
  if (ev.GetEndShower()) {
    return;
  }

  // Set to False, only set to True if accpeted
  ev.SetAcceptEmission(false);

  curandState state = states[idx];

  double z = 0.0;
  double zp = ev.GetWinParam(0);
  double zm = 1.0 - zp;

  int kernel = ev.GetWinSF();
  int ij_flav = ev.GetParton(ev.GetWinDipole(0)).GetPid();

  double rand = curand_uniform(&state);
  states[idx] = state;

  // Generate z
  if (kernel < 10 && ij_flav != 21) {
    z = 1.0 + (zp - 1.0) * pow((1.0 - zm) / (1.0 - zp), rand);
  } else if (kernel < 11) {
    z = 1.0 + (zp - 1.0) * pow((1.0 - zm) / (1.0 - zp), rand);
  } else if (kernel < 16) {
    z = zm + (zp - zm) * rand;
  }

  double y = ev.GetShowerT() / ev.GetWinParam(1) / z / (1.0 - z);

  double f = 0.0;
  double g = 0.0;
  double value = 0.0;
  double estimate = 0.0;

  // CS Kernel: y can't be 1
  if (y < 1.0) {
    // Value and Estimate
    if (kernel < 10 && ij_flav != 21) {
      value = kCF * (2.0 / (1.0 - z * (1.0 - y)) - (1.0 + z));
      estimate = kCF * 2.0 / (1.0 - z);
    } else if (kernel < 11) {
      value = kCA / 2.0 * (2.0 / (1.0 - z * (1.0 - y)) - 2.0 + z * (1.0 - z));
      estimate = kCA / (1.0 - z);
    } else if (kernel < 16) {
      value = kTR / 2.0 * (1.0 - 2.0 * z * (1.0 - z));
      estimate = kTR / 2.0;
    }

    f = (1.0 - y) * ev.GetAsVeto() * value;
    g = asmax * estimate;

    if (curand_uniform(&state) < f / g) {
      ev.SetAcceptEmission(true);
      ev.SetShowerZ(z);
      ev.SetShowerY(y);
    }
    states[idx] = state;
  }
}

// -----------------------------------------------------------------------------

// Do Splitting
__global__ void doSplitting(Event *events, curandState *states, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event &ev = events[idx];

  // Do not run if the shower has ended
  if (ev.GetEndShower()) {
    return;
  }

  if (!ev.GetAcceptEmission()) {
    return;
  }

  curandState state = states[idx];

  double phi = 2.0 * M_PI * curand_uniform(&state);
  states[idx] = state;

  int win_ij = ev.GetWinDipole(0);
  int win_k = ev.GetWinDipole(1);

  // Make Kinematics
  Vec4 moms[3] = {Vec4(), Vec4(), Vec4()};

  MakeKinematics(moms, ev.GetShowerZ(), ev.GetShowerY(), phi,
                 ev.GetParton(win_ij).GetMom(), ev.GetParton(win_k).GetMom());

  // Adjust Colors

  // ---------------------------------------------------------------------------
  /**
   * I think there is a more elegant way of handling this, like a __device__
   * function, but I'll leave it like this for now (Feb 2024), until paper
   * is done.
   */

  // Get Flavs from Kernel Number
  int kernel = ev.GetWinSF();

  // flavour: -5 -4 -3 -2 -1  1  2  3  4  5
  // index:    0  1  2  3  4  5  6  7  8  9
  int flavs[3];
  if (kernel < 10 && ev.GetParton(ev.GetWinDipole(0)).GetPid() != 21) {
    if (kernel < 5) {
      flavs[0] = kernel - 5;
      flavs[1] = kernel - 5;
      flavs[2] = 21;
    } else if (kernel < 10) {
      flavs[0] = kernel - 4;
      flavs[1] = kernel - 4;
      flavs[2] = 21;
    }
  } else if (kernel < 11) {
    flavs[0] = 21;
    flavs[1] = 21;
    flavs[2] = 21;
  } else if (kernel < 16) {
    flavs[0] = 21;
    flavs[1] = kernel - 10;
    flavs[2] = -1 * (kernel - 10);
  }

  // ---------------------------------------------------------------------------

  int colij[2] = {ev.GetParton(win_ij).GetCol(),
                  ev.GetParton(win_ij).GetAntiCol()};

  int colk[2] = {ev.GetParton(win_k).GetCol(),
                 ev.GetParton(win_k).GetAntiCol()};

  int coli[2] = {0, 0};
  int colj[2] = {0, 0};

  double rand = curand_uniform(&state);
  states[idx] = state;

  MakeColours(ev, coli, colj, flavs, colij, colk, rand);

  // Modify Splitter
  ev.SetPartonPid(win_ij, flavs[1]);
  ev.SetPartonMom(win_ij, moms[0]);
  ev.SetPartonCol(win_ij, coli[0]);
  ev.SetPartonAntiCol(win_ij, coli[1]);

  // Modify Recoiled Spectator
  ev.SetPartonMom(win_k, moms[2]);

  // Add Emitted Parton
  Parton em = Parton(flavs[2], moms[1], colj[0], colj[1]);
  ev.SetParton(ev.GetSize(), em);

  // Increment Emissions (IMPORTANT)
  ev.IncrementEmissions();
}

// -----------------------------------------------------------------------------

/*
__global__ void countBools(Event *events, int *trueCount,
                           int *falseCount, int N) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= N) {
    return;
  }

  Event &ev = events[idx];

  if (ev.GetEndShower()) {
    return;
  }

  if (!ev.GetAcceptEmission()){
    atomicAdd(trueCount, 1);
  } else {
    atomicAdd(falseCount, 1);
  }
}
*/

__device__ void MakeKinematics(Vec4 *kinematics, const double z, const double y,
                               const double phi, const Vec4 pijt,
                               const Vec4 pkt) {
  Vec4 Q = pijt + pkt;

  // Generating the Momentum (0, kt1, kt2, 0)
  double rkt = sqrt(Q.M2() * y * z * (1. - z));

  Vec4 kt1 = pijt.Cross(pkt);
  if (kt1.P() < 1.e-6) {
    Vec4 xaxis(0., 1., 0., 0.);
    kt1 = pijt.Cross(xaxis);
  }
  kt1 = kt1 * (rkt * cos(phi) / kt1.P());

  Vec4 kt2cms = Q.Boost(pijt);
  kt2cms = kt2cms.Cross(kt1);
  kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.P());
  Vec4 kt2 = Q.BoostBack(kt2cms);

  // Conversion to {i, j, k} basis
  Vec4 pi = pijt * z + pkt * ((1. - z) * y) + kt1 + kt2;
  Vec4 pj = pijt * (1. - z) + pkt * (z * y) - kt1 - kt2;
  Vec4 pk = pkt * (1. - y);

  // No need to do *kinematics[0], for arrays the elements are already
  // pointers
  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
}

// -----------------------------------------------------------------------------

void runShower(thrust::device_vector<Event> &d_events) {
  // Number of Events - Can get from d_events.size()
  int N = d_events.size();

  // Set up the device alphaS
  AlphaS *d_as;
  cudaMalloc(&d_as, sizeof(AlphaS));
  asSetupKernel<<<1, 1>>>(d_as, mz, asmz);
  syncGPUAndCheck("asSetupKernel");

  // Allocate device memory for completed events counter
  int *d_completed;
  cudaMalloc(&d_completed, sizeof(int));
  cudaMemset(d_completed, 0, sizeof(int));

  // Allocate space for curand states
  curandState *d_states;
  cudaMalloc(&d_states, N * sizeof(curandState));

  // Initialize the states
  initCurandStates<<<(N + 255) / 256, 256>>>(d_states, N);

  // Store the number of finished events per cycle
  std::vector<int> completedPerCycle;

  // Use a pointer to the device events
  Event *d_events_ptr = thrust::raw_pointer_cast(d_events.data());

  // ---------------------------------------------------------------------------
  // Prepare the Shower

  DEBUG_MSG("Running @prepShower");
  prepShower<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("prepShower");

  // ---------------------------------------------------------------------------
  // Run the Shower

  int completed = 0;
  int cycle = 0;
  while (completed < N) {
    // Run all the kernels here...
    // -------------------------------------------------------------------------
    // Select the winner kernel

    DEBUG_MSG("Running @selectWinnerSplitFunc");
    selectWinnerSplitFunc<<<(N + 255) / 256, 256>>>(d_events_ptr, d_states, N);
    syncGPUAndCheck("selectWinnerSplitFunc");

    // -------------------------------------------------------------------------
    // Check Cutoff

    DEBUG_MSG("Running @checkCutoff");
    checkCutoff<<<(N + 255) / 256, 256>>>(d_events_ptr, d_completed, tC, N);
    syncGPUAndCheck("checkCutoff");

    // -------------------------------------------------------------------------
    // Calculate AlphaS for Veto Algorithm

    DEBUG_MSG("Running @asShowerKernel");
    asShowerKernel<<<(N + 255) / 256, 256>>>(d_as, d_events_ptr, N);
    syncGPUAndCheck("asShowerKernel");

    // -------------------------------------------------------------------------
    // Veto Algorithm

    DEBUG_MSG("Running @vetoAlg");
    vetoAlg<<<(N + 255) / 256, 256>>>(d_events_ptr, d_states, N);
    syncGPUAndCheck("vetoAlg");

    // -------------------------------------------------------------------------
    // Splitting Algorithm

    DEBUG_MSG("Running @doSplitting");
    doSplitting<<<(N + 255) / 256, 256>>>(d_events_ptr, d_states, N);
    syncGPUAndCheck("doSplitting");

    // -------------------------------------------------------------------------
    // Import the Number of Completed Events

    cudaMemcpy(&completed, d_completed, sizeof(int), cudaMemcpyDeviceToHost);
    cycle++;

    // Until Paper is Published, we will use this
    completedPerCycle.push_back(completed);

    // -------------------------------------------------------------------------
    // Print number of Accepted / Vetoed Events - for A. V.

    /*
    // TRUE means that the event is vetoed
    // FALSE means that the event is accepted

    int *d_trueCount, *d_falseCount;
    cudaMalloc(&d_trueCount, sizeof(int));
    cudaMalloc(&d_falseCount, sizeof(int));
    cudaMemset(d_trueCount, 0, sizeof(int));
    cudaMemset(d_falseCount, 0, sizeof(int));

    DEBUG_MSG("Running @countBools");
    countBools<<<(N + 255) / 256, 256>>>(d_events_ptr, d_trueCount,
                                         d_falseCount, N);
    syncGPUAndCheck("countBools");

    int h_trueCount(0), h_falseCount(0);  // Number of vetoed events
    cudaMemcpy(&h_trueCount, d_trueCount, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_falseCount, d_falseCount, sizeof(int),
               cudaMemcpyDeviceToHost);

    std::cout << cycle << ", " << N - completed << ", " << h_trueCount << ", "
              << h_falseCount << std::endl;
    */
  }

  // ---------------------------------------------------------------------------
  // Write completedPerCycle to file
  std::ofstream file("gaps-cycles.dat");
  for (auto &i : completedPerCycle) {
    file << i << std::endl;
  }

  // ---------------------------------------------------------------------------
  // Clean Up Device Memory
  cudaFree(d_completed);
}
