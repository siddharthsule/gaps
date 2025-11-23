#include "mczinc.cuh"

__global__ void calculate_mczinc(const event* events, double* results, int n) {
  /**
   * @brief Calculate the Z Boson Variables, akin to rivet's MC_ZINC
   *
   * @param ev The event object
   * @param results The array to store the results
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Observables Preamble
  if (!events[idx].get_validity()) return;
  const event& ev = events[idx];
  // ---------------------------------------------

  // calculate the z boson momentum
  // LO -> e+ e-, NLO -> Z
  vec4 z_mom;
  if (ev.get_particle(2).get_pid() == 23) {
    z_mom = ev.get_particle(2).get_mom();
  } else {
    z_mom = ev.get_particle(2).get_mom() + ev.get_particle(3).get_mom();
  }

  // set the z boson observables (m, pt, phi, y)
  results[20 * idx + 0] = z_mom.m();
  results[20 * idx + 1] = z_mom.pt();
  results[20 * idx + 2] = z_mom.phi();
  results[20 * idx + 3] = z_mom.rapidity();

  // set the lepton observables (pt, y)
  if (!(ev.get_particle(2).get_pid() == 23)) {
    results[20 * idx + 4] = ev.get_particle(2).get_mom().pt();
    results[20 * idx + 5] = ev.get_particle(3).get_mom().pt();
    results[20 * idx + 6] = ev.get_particle(2).get_mom().eta();
    results[20 * idx + 7] = ev.get_particle(3).get_mom().eta();
  }
}