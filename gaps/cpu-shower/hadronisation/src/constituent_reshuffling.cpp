#include "hadronisation.h"

double hadronisation::f_reshuffling(double k, double* masses, double ecms,
                                    event& ev) const {
  /**
   * @brief function to find the factor k such that the sum of the reshuffled
   * momenta equals Ecms
   *
   * for i in range(2, len(event)):
   *     total += m.sqrt(k**2 * event[i].mom.P2() + masses[i-2]**2)
   *
   * @param k the factor to find
   * @param masses the constituent masses of the particles
   * @param ev the event to reshuffle
   * @return the difference between the sum of the reshuffled momenta and Ecms
   */

  // Total
  double total = 0.;

  // Loop over final state particles
  for (int i = 2; i < ev.get_size(); i++) {
    // Get the three momentum squared of the particle
    double p2 = ev.get_particle(i).get_mom().p2();
    double m2 = masses[i - 2] * masses[i - 2];

    // Add to the total
    total += sqrt(k * k * p2 + m2);
  }
  return total - ecms;
}

void hadronisation::constituent_reshuffling(event& ev) const {
  /**
   * @brief reshuffle the quarks and gluons to their constituent masses
   *
   * Reshuffles the momenta of the massless particles from the parton shower to
   * massive particles. Fortunately, the sum of the momenta of the massless
   * particles is (Ecms, 0), and we don't need to worry about boosting. We
   * mainly need to find a factor k such that:
   *
   * Sum( sqrt(k^2*|p|_i^2 + m_i^2) ) = Ecms
   *
   * where p_i is the momentum of the massless particle, |p|_i is the
   * three momentum of the massless particle, and m_i is the constituent mass of
   * the massive particle. then, we can set the momentum of the massive
   * particle to:
   *
   * k0 = sqrt(k^2*|p|_i^2 + m_i^2), |k|_i = |p|_i * k
   *
   * @param ev the event to reshuffle
   */

  // Calculate total momentum sum of all particles
  vec4 total_momentum;

  // Get a list of constituent masses (length = event size - 2 initial)
  double masses[max_particles];

  // Skip the first two particles (incoming beams)
  for (int i = 2; i < ev.get_size(); i++) {
    // Add to the total momentum
    total_momentum = total_momentum + ev.get_particle(i).get_mom();

    // Append the correct constituent mass
    if (ev.get_particle(i).get_pid() == 21) {
      masses[i - 2] = gluon_mass;  // gluon
    } else {
      int fl = abs(ev.get_particle(i).get_pid());  // quark
      masses[i - 2] = const_mass[fl - 1];
    }
  }

  // Define the total energy as the center-of-mass energy (Ecms)
  double ecms = total_momentum[0];

  // Do Bisection method to find k
  double klo = 0.0;
  double khi = 1.0;
  for (int iter = 0; iter < 50; iter++) {
    double kmid = 0.5 * (klo + khi);
    if (f_reshuffling(kmid, masses, ecms, ev) > 0.0) {
      khi = kmid;
    } else {
      klo = kmid;
    }
    if (khi - klo < 1e-12) {
      break;
    }
  }
  double k = 0.5 * (klo + khi);

  // Reshuffle the momentum of the particles
  for (int i = 2; i < ev.get_size(); i++) {
    // Get the three momentum squared of the particle
    double p2 = ev.get_particle(i).get_mom().p2();
    double m2 = masses[i - 2] * masses[i - 2];

    // Reshuffle the momentum of the particle
    double e = sqrt(k * k * p2 + m2);
    double px = k * ev.get_particle(i).get_mom()[1];
    double py = k * ev.get_particle(i).get_mom()[2];
    double pz = k * ev.get_particle(i).get_mom()[3];

    // Set the new momentum of the particle
    ev.set_particle_mom(i, vec4(e, px, py, pz));
  }
}
