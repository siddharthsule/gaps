#include "mczinc.h"

void calculate_mczinc(const event& ev, double* results) {
  /**
   * @brief Calculate the Z Boson Variables, akin to rivet's MC_ZINC
   *
   * @param ev The event object
   * @param results The array to store the results
   */

  if (!ev.get_validity()) {
    return;
  }

  // calculate the z boson momentum
  // LO -> e+ e-, NLO -> Z
  vec4 z_mom;
  if (ev.get_particle(2).get_pid() == 23) {
    z_mom = ev.get_particle(2).get_mom();
  } else {
    z_mom = ev.get_particle(2).get_mom() + ev.get_particle(3).get_mom();
  }

  // set the z boson observables (m, pt, phi, y)
  results[0] = z_mom.m();
  results[1] = z_mom.pt();
  results[2] = z_mom.phi();
  results[3] = z_mom.rapidity();

  // set the lepton observables (pt, y)
  if (!(ev.get_particle(2).get_pid() == 23)) {
    results[4] = ev.get_particle(2).get_mom().pt();
    results[5] = ev.get_particle(3).get_mom().pt();
    results[6] = ev.get_particle(2).get_mom().eta();
    results[7] = ev.get_particle(3).get_mom().eta();
  }
}