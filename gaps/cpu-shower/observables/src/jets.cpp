#include "jets.h"

// -----------------------------------------------------------------------------
// Durham Algorithm for LEP

double yij(const vec4& p, const vec4& q, double ecm2) {
  /**
   * @brief Calculate the yij value for the Durham Algorithm
   *
   * @param p The 4-momentum of the first particle
   * @param q The 4-momentum of the second particle
   * @param ecm2 The center of mass energy squared
   * @return double The yij value
   */

  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = std::min(p[0], q[0]);
  double max_pq = std::max(pq / std::sqrt(p.p2() * q.p2()), -1.);
  return 2. * std::pow(min_pq, 2) * (1. - std::min(max_pq, 1.)) / ecm2;
}

void cluster_durham(const event& ev, double* results) {
  /**
   * @brief Cluster the event int jets using the Durham algorithm for LEP
   *
   * @param ev The event object
   * @param results The array to store the results
   */

  if (!ev.get_validity()) {
    return;
  }

  // get the center of mass energy squared
  double ecm2 =
      (ev.get_particle(0).get_mom() + ev.get_particle(1).get_mom()).m2();

  // extract the 4-momenta of the particles
  vec4 p[max_particles];
  for (int i = 2; i < ev.get_size(); ++i) {
    p[i - 2] = ev.get_particle(i).get_mom();
  }

  // kt2 will store the kt2 values for each clustering step
  // if not changed, set to -1 so we can ignore when histogramming
  double kt2[max_particles] = {-1.};
  int counter = 0;

  // num particles (which will change when clustered)
  int n_particles = ev.get_size() - 2;

  // imap will store the indices of the particles
  int imap[max_particles];
  for (int i = 0; i < ev.get_size() - 2; ++i) {
    imap[i] = i;
  }

  // kt2ij will store the kt2 values for each pair of particles
  double kt2ij[max_particles][max_particles] = {0.};
  double dmin = 1.;
  int ii = 0, jj = 0;
  for (int i = 0; i < n_particles; ++i) {
    for (int j = 0; j < i; ++j) {
      double dij = kt2ij[i][j] = yij(p[i], p[j], ecm2);

      // Find the smallest Yij measure
      if (dij < dmin) {
        dmin = dij;
        ii = i;
        jj = j;
      }
    }
  }

  // cluster the particles
  // n_particles = 2 as the event is a 2-jet event, and we are only interested
  // in y23...
  while (n_particles > 2) {
    // Our observable is the dmin required to combine the particles
    // with the smallest Yij measure. For example, y23 is the smallest
    // Yij measure required to combine a 3 particle event into a 2
    // jet event.
    kt2[counter] = dmin;

    // Let's add an example to illustrate the following steps:
    // Eg: [0, 1, 2, 3, 4] with the smallest Yij between 1 and 2
    // ii = 1, jj = 2, imap = [0, 1, 2, 3, 4]

    // jjx = index of the second particle in the pair (ii, jj)
    // Eg: jjx = jj = 2
    int jjx = imap[jj];

    // The momentum of ii is added to jjx
    // Eg: p[2] += p[1]
    p[jjx] = p[jjx] + p[imap[ii]];

    // Reduce the number of particles by 1
    --n_particles;
    counter++;

    // Remove the particle with the index ii by shifting the indices
    // Eg : [0, 1, 2, 3, 4]->[0, 2, 3, 4, 4]
    for (int i = ii; i < n_particles; ++i) {
      imap[i] = imap[i + 1];
    }

    // Update the Yij measure matrix. Instead of recalculating all the
    // Yij measures, we only need to update the Yij measures of the
    // particles that were combined with the particle jjx

    // Note that using imap, we skip the particle with the index ii, as
    // it was "removed" from the event

    // Eg: jj = 2, jjx = 2, imap = [0, 2, 3, 4, 4]
    // for j in range(2): kt2ij[2][imap[j]] = self.Yij(p[2], p[imap[j]])
    // kt2ij[2][imap[0]] = Yij(p[2], p[imap[0]]) -> kt2ij[2][0]  = Yij(p[2],
    // p[0]) kt2ij[2][imap[1]] = Yij(p[2], p[imap[1]]) -> kt2ij[2][2]  =
    // Yij(p[2], p[2])
    for (int j = 0; j < jj; ++j) {
      kt2ij[jjx][imap[j]] = yij(p[jjx], p[imap[j]], ecm2);
    }

    // Eg : jj = 2, jjx = 2, imap = [0, 2, 3, 4, 4]
    // for i in range(3, 4) : kt2ij[imap[i]][2] = self.Yij(p[2], p[imap[i]])
    // kt2ij[3][2] = Yij(p[2], p[3])
    // kt2ij[4][2] = Yij(p[2], p[4])
    for (int i = jj + 1; i < n_particles; ++i) {
      kt2ij[imap[i]][jjx] = yij(p[jjx], p[imap[i]], ecm2);
    }

    // Find the next smallest Yij measure
    dmin = 1.;
    for (int i = 0; i < n_particles; ++i) {
      for (int j = 0; j < i; ++j) {
        // Get the Yij measure.Because of the updated imap, we
        // skip the particle with the index ii in this step
        double dij = kt2ij[imap[i]][imap[j]];

        // Find the smallest Yij measure
        if (dij < dmin) {
          dmin = dij;
          ii = i;
          jj = j;
        }
      }
    }
  }

  // store the results (y23, y34, y45, y56, tvalue, tzoomd, hjm, ljm, wjb, njb)
  results[0] = counter > 0 ? log10(kt2[counter - 1 - 0]) : -50.;
  results[1] = counter > 1 ? log10(kt2[counter - 1 - 1]) : -50.;
  results[2] = counter > 2 ? log10(kt2[counter - 1 - 2]) : -50.;
  results[3] = counter > 3 ? log10(kt2[counter - 1 - 3]) : -50.;

  return;
}