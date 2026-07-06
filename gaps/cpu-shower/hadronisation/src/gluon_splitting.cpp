#include "hadronisation.h"

void hadronisation::force_gluons_to_split(event& ev) const {
  /**
   * @brief Force gluons to split isotropically into quark-antiquark pairs
   *
   * The allowed flavours of the q-qbar pair depend on whether the mass of
   * the gluon is twice the mass of the quark. Usually, this means only
   * up down and strange quarks can be produced from gluons, while charm
   * and bottom quarks can only be produced from gluons with sufficient
   * energy.
   */

  // Loop over final state particles in the event
  // Cache size - don't process newly added antiquarks
  int n_particles = ev.get_size();
  for (int i = 2; i < n_particles; i++) {
    // If Gluon
    if (ev.get_particle(i).get_pid() != 21) continue;

    // Get the colour indices of the gluon
    int col = ev.get_particle(i).get_col();
    int acol = ev.get_particle(i).get_acol();

    // Check for colour-singlet gluon
    if (col == acol) {
      std::cerr << "gluon_splitting: gluon has same colour indices!"
                << std::endl;
      continue;
    }

    // Pick flavour, u and d more likely than s
    int fl = select_qq_flavour(ev.gen_random());

    // Define masses of the quarks
    double m = const_mass[fl - 1];

    // Generate Kinematics for the gluon splitting
    vec4 p1, p2;
    double rho_1 = ev.gen_random();
    double rho_2 = ev.gen_random();
    kallen(ev.get_particle(i).get_mom(), m, m, p1, p2, rho_1, rho_2);

    // Set the current particle as the quark
    ev.set_particle_pid(i, fl);
    ev.set_particle_mom(i, p1);
    ev.set_particle_col(i, col);
    ev.set_particle_acol(i, 0);

    // Add the antiquark with opposite colour
    ev.add_emission(particle(-fl, p2, 0, acol));
  }

  // Check if there are any gluons left
  for (int i = 2; i < ev.get_size(); i++) {
    if (ev.get_particle(i).get_pid() == 21)
      std::cerr << "gluon_splitting: gluon not split!" << std::endl;
  }
}