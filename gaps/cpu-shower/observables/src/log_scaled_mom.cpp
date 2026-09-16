#include "log_scaled_mom.h"

// particle_to_charge, the same charged test the multiplicity applies
#include "chargedmult.h"

void fill_log_scaled_mom(const event& ev, histo1d& hist, double weight) {
  /**
   * @brief Fill the Log Scaled Momentum of every charged particle for LEP
   *
   * @param ev The event object
   * @param hist The histogram to fill
   * @param weight The event weight
   */

  if (!ev.get_validity()) {
    return;
  }

  // Get the beam momentum
  double p_beam = ev.get_particle(0).get_mom().p();
  if (p_beam <= 0.) {
    return;
  }

  // Compute the scaled momentum
  for (int i = 2; i < ev.get_size(); ++i) {
    // Charged only
    if (particle_to_charge(ev.get_particle(i).get_pid()) == 0) {
      continue;
    }

    // A particle at rest has no log scaled momentum
    double p = ev.get_particle(i).get_mom().p();
    if (p <= 0.) {
      continue;
    }

    // Add to bin
    hist.fill(-log(p / p_beam), weight);
  }
}
