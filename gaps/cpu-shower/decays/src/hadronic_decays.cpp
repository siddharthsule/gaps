#include "hadronic_decays.h"

// -----------------------------------------------------------------------------
// Stability Check

bool hadronic_decays::is_stable(int pid) const {
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

bool hadronic_decays::decay_one(event& ev, int i) const {
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
    std::cerr << "hadronic_decays: " << parent.get_pid() << " of mass "
              << parent.get_mom().m() << " GeV left undecayed" << std::endl;
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

void hadronic_decays::run(event& ev) const {
  /**
   * @brief Decay every unstable particle until the event is final.
   *
   * @param ev the event to decay
   */
  // Check for overflow
  if (ev.get_overflowed()) return;

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

  std::cerr << "hadronic_decays: event still holds unstable particles after "
            << max_generations << " generations" << std::endl;
}
