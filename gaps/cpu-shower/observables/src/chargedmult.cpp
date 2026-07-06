#include "chargedmult.h"

#include "hadrons.h"

void calculate_chargedmult(const event& ev, double* result) {
  /**
   * @brief Calculate the Charged Multiplicity for LEP
   *
   * @param ev The event object
   * @param result The array to store the result
   */

  if (!ev.get_validity()) {
    return;
  }

  int nch = 0;
  for (int i = 2; i < ev.get_size(); ++i) {
    if (particle_to_charge(ev.get_particle(i).get_pid()) != 0) {
      ++nch;
    }
  }

  result[0] = static_cast<double>(nch);
}