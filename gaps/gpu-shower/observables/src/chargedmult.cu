#include "chargedmult.cuh"

__global__ void calculate_chargedmult(const event* events, double* results,
                                      int result_slot, int n) {
  /**
   * @brief Calculate the Charged Multiplicity for LEP
   *
   * @param events The events array
   * @param results The array to store the results
   * @param result_slot The slot in the results array to store the result
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

  int nch = 0;
  for (int i = 2; i < ev.get_size(); ++i) {
    if (particle_to_charge(ev.get_particle(i).get_pid()) != 0) ++nch;
  }

  results[20 * idx + result_slot] = static_cast<double>(nch);
}
