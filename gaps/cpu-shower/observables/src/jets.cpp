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
  double min_pq = min(p[0], q[0]);
  double max_pq = max(pq / sqrt(p.p2() * q.p2()), -1.);
  return 2. * pow(min_pq, 2) * (1. - min(max_pq, 1.)) / ecm2;
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

// -----------------------------------------------------------------------------
// Generalized kt Algorithm for LHC

void bubble_sort_pt(vec4* moms, int n_moms) {
  /**
   * @brief Bubble Sort the Momenta in Descending Order of pt()
   *
   * @param moms The array of momenta
   * @param n_moms The number of momenta
   */

  if (n_moms == 0) return;

  for (int i = 0; i < n_moms - 1; i++) {
    for (int j = 0; j < n_moms - i - 1; j++) {
      if (moms[j].pt() < moms[j + 1].pt()) {
        vec4 temp = moms[j];
        moms[j] = moms[j + 1];
        moms[j + 1] = temp;
      }
    }
  }
}

double dR2(const vec4& p, const vec4& q) {
  /**
   * @brief Calculate the dR^2 value for the Gen-kt Algorithm
   *
   * @param p The 4-momentum of the first particle
   * @param q The 4-momentum of the second particle
   * @return double The dR^2 value
   */

  // Delta R_{ij}^2 = (y_i - y_j)^2 + (phi_i - phi_j)^2
  // double dy = p.rapidity() - q.rapidity();
  double dy = p.eta() - q.eta();
  double dphi = p.delta_phi(q);
  return dy * dy + dphi * dphi;
}

double dij(const vec4& p, const vec4& q) {
  /**
   * @brief Calculate the dij value for the Gen-kt Algorithm
   *
   * @param p The 4-momentum of the first particle
   * @param q The 4-momentum of the second particle
   * @param power The power of the distance measure
   * @param R The jet parameter
   * @return double The dij value
   */

  // p_prime = pt^(2*power) =(pt^2)^power
  double p_prime = pow(p.pt2(), power);
  double q_prime = pow(q.pt2(), power);

  // dij measure
  return min(p_prime, q_prime) * dR2(p, q) / (R * R);
}

void cluster_genkt(const event& ev, double* results) {
  /**
   * @brief Cluster the event int jets using the Gen-kt algorithm for LHC
   *
   * @param ev The event object
   * @param results The array to store the results
   */

  if (!ev.get_validity()) {
    return;
  }

  // The particle list is either:
  // LO: p(q), p(qbar), e+, e-
  // NLO: p(q), p(qbar), Z
  int start = ev.get_particle(2).get_pid() == 23 ? 3 : 4;

  // extract the 4-momenta of the particles
  vec4 p[max_particles];
  for (int i = start; i < ev.get_size(); ++i) {
    p[i - start] = ev.get_particle(i).get_mom();
  }

  // num particles (which will change when clustered)
  int n_particles = ev.get_size() - start;

  // Array to store jets: vec4 objects (say 10 for safety?)
  const int max_jets = 10;
  vec4 jets[max_jets];
  int jet_counter = 0;

  // imap will store the indices of the particles
  int imap[max_particles];
  for (int i = 0; i < ev.get_size() - start; ++i) {
    imap[i] = i;
  }

  // dmax is an arbitrary large number
  double dmax = 1e6;

  // dij_matrix will store the dij values for each pair of particles
  double dij_matrix[max_particles][max_particles] = {-1.};
  double diB[max_particles] = {-1.};

  // Define a variable to store whether dmin comes from dij_matrix or diB
  bool dmin_is_diB = false;

  // Calculate the dij values for each pair of particles
  double dmin = dmax;
  int ii = 0, jj = 0;
  for (int i = 0; i < n_particles; ++i) {
    // Calculate the dij measure
    for (int j = 0; j < i; ++j) {
      dij_matrix[i][j] = dij(p[i], p[j]);

      // Check if dij is the smallest
      if (dij_matrix[i][j] < dmin) {
        dmin = dij_matrix[i][j];
        ii = i;
        jj = j;
        dmin_is_diB = false;
      }
    }

    // Calculate the diB measure
    diB[i] = pow(p[i].pt2(), power);

    // Check if diB is the smallest
    if (diB[i] < dmin) {
      dmin = diB[i];
      ii = i;
      dmin_is_diB = true;
    }
  }

  // cluster the particles
  while (n_particles > 0) {
    // Break when all particles are combined and jets are stored
    if ((dmin == dmax) || (jet_counter == max_jets)) {
      break;
    }

    // printf("Loop Iteration:\n");
    // printf("dmin = %f\n", dmin);
    // printf("ii = %d\n", ii);
    // printf("jj = %d\n", jj);
    // if (dmin_is_diB) {
    //   printf("dmin is from diB\n");
    // } else {
    //   printf("dmin is from dij_matrix\n");
    // }
    // printf("n_particles = %d\n", n_particles);
    // printf("jet_counter = %d\n", jet_counter);

    // If dmin == diB, then we have a final jet. Store the jet and
    // remove the particle from the event
    if (dmin_is_diB) {
      // Let's add an example to illustrate the following steps:
      // Eg: [0, 1, 2, 3, 4] with the dmin coming from diB of 1
      // ii = 1, imap = [0, 1, 2, 3, 4]

      // Store the jet
      jets[jet_counter] = p[imap[ii]];
      jet_counter++;

      // Reduce the number of particles by 1
      --n_particles;

      // Remove the particle with the index ii by shifting the indices
      // Eg: [0, 1, 2, 3, 4] -> [0, 2, 3, 4, 4]
      for (int i = ii; i < n_particles; ++i) {
        imap[i] = imap[i + 1];
      }

      // The particle is removed from imap, so we don't need to update
      // dij or diB - this element will never be accessed again!
    }

    else {
      // Combine the two particles with the smallest dij measure
      // Let's add an example to illustrate the following steps:
      // Eg: [0, 1, 2, 3] with the smallest dij between 1 and 2
      // ii = 1, jj = 2, imap = [0, 1, 2, 3]

      // jjx = index of the second particle in the pair (ii, jj)
      // Eg: jjx = jj = 2
      int jjx = imap[jj];

      // The momentum of ii is added to jjx
      // Eg: p[2] += p[1]
      p[jjx] = p[jjx] + p[imap[ii]];

      // Reduce the number of particles by 1
      --n_particles;

      // Remove the particle with the index ii by shifting the indices
      // Eg : [0, 1, 2, 3, 4]->[0, 2, 3, 4, 4]
      for (int i = ii; i < n_particles; ++i) {
        imap[i] = imap[i + 1];
      }

      // Update the dij measure matrix. Instead of recalculating all the
      // dij measures, we only need to update the dij measures of the
      // particles that were combined with the particle jjx

      // Note that using imap, we skip the particle with the index ii, as
      // it was "removed" from the event

      // Eg: jj = 2, jjx = 2, imap = [0, 2, 3, 4, 4]
      // for j in range(2): dij[2][imap[j]] = self.dij(p[2], p[imap[j]])
      // dij[2][imap[0]] = dij(p[2], p[imap[0]]) -> dij[2][0]  = dij(p[2],p[0])
      // dij[2][imap[1]] = dij(p[2], p[imap[1]]) -> dij[2][2]  = dij(p[2],p[2])
      for (int j = 0; j < jj; ++j) {
        dij_matrix[jjx][imap[j]] = dij(p[jjx], p[imap[j]]);
      }

      // Eg : jj = 2, jjx = 2, imap = [0, 2, 3, 4, 4]
      // for i in range(3, 4) : dij[imap[i]][2] = self.dij(p[2], p[imap[i]])
      // dij[3][2] = dij(p[2], p[3])
      // dij[4][2] = dij(p[2], p[4])
      for (int i = jj + 1; i < n_particles; ++i) {
        dij_matrix[imap[i]][jjx] = dij(p[jjx], p[imap[i]]);
      }

      // Update the PT2 element of the particle jjx
      // Eg: jj = 2, jjx = 2, imap = [0, 2, 3, 3]
      // diB[2] = pt^(2*power)(p[2])
      diB[jjx] = pow(p[jjx].pt2(), power);
    }

    // Find the next smallest dij measure
    dmin = dmax;
    for (int i = 0; i < n_particles; ++i) {
      for (int j = 0; j < i; ++j) {
        // Get the dij measure.Because of the updated imap, we
        // skip the particle with the index ii in this step
        double dij = dij_matrix[imap[i]][imap[j]];

        // Find the smallest dij measure
        if (dij < dmin) {
          dmin = dij;
          ii = i;
          jj = j;
          dmin_is_diB = false;
        }
      }

      // Check if the existing diB is the smallest
      if (diB[imap[i]] < dmin) {
        dmin = diB[imap[i]];
        ii = i;
        dmin_is_diB = true;
      }
    }
  }

  // Sort the jets by PT
  bubble_sort_pt(jets, jet_counter);

  // Optionally, sort the particles by pt and call them jets
  // n_particles = ev.get_size() - 4;
  // bubble_sort_pt(p, n_particles);
  // for (int i = 0; i < max_jets; ++i) {
  //   jets[i] = p[i];
  // }
  // jet_counter = n_particles < max_jets ? n_particles : max_jets;

  // Remove jets with PT < 0 GeV
  while (jet_counter > 0 && jets[jet_counter - 1].pt() <= 0.) {
    jet_counter--;
  }

  // Calculate the z boson momentum
  // LO -> e+ e-, NLO -> Z
  vec4 z_mom;
  if (ev.get_particle(2).get_pid() == 23) {
    z_mom = ev.get_particle(2).get_mom();
  } else {
    z_mom = ev.get_particle(2).get_mom() + ev.get_particle(3).get_mom();
  }

  double z_j1_eta, z_j1_dR;
  if (jet_counter > 0) {
    z_j1_eta = z_mom.eta() - jets[0].eta();
    z_j1_dR = sqrt(dR2(z_mom, jets[0]));
  } else {
    z_j1_eta = -50.;
    z_j1_dR = -50.;
  }

  double j1_j2_dR, j1_j3_dR, j2_j3_dR;
  if (jet_counter > 1) {
    j1_j2_dR = sqrt(dR2(jets[0], jets[1]));
  } else {
    j1_j2_dR = -50.;
  }

  if (jet_counter > 2) {
    j1_j3_dR = sqrt(dR2(jets[0], jets[2]));
    j2_j3_dR = sqrt(dR2(jets[1], jets[2]));
  } else {
    j1_j3_dR = -50.;
    j2_j3_dR = -50.;
  }

  // store the first three jets' pt and eta
  double jet1_pt = jet_counter > 0 ? jets[0].pt() : -50.;
  double jet1_eta = jet_counter > 0 ? jets[0].eta() : -50.;
  double jet2_pt = jet_counter > 1 ? jets[1].pt() : -50.;
  double jet2_eta = jet_counter > 1 ? jets[1].eta() : -50.;
  double jet3_pt = jet_counter > 2 ? jets[2].pt() : -50.;
  double jet3_eta = jet_counter > 2 ? jets[2].eta() : -50.;

  // Jet PT (with eta cuts)
  results[8] = (jet_counter > 0 && fabs(jet1_eta) < 5.0) ? jet1_pt : -50.;
  results[9] = (jet_counter > 1 && fabs(jet2_eta) < 5.0) ? jet2_pt : -50.;
  results[10] = (jet_counter > 2 && fabs(jet3_eta) < 5.0) ? jet3_pt : -50.;

  // Jet Eta (with pt cuts)
  results[11] = (jet_counter > 0 && jet1_pt > 5.0) ? jet1_eta : -50.;
  results[12] = (jet_counter > 1 && jet2_pt > 5.0) ? jet2_eta : -50.;
  results[13] = (jet_counter > 2 && jet3_pt > 5.0) ? jet3_eta : -50.;

  // Count jets satisfying |eta| < 5 and pt > 5
  int n_good_jets = 0;
  for (int i = 0; i < jet_counter; ++i) {
    if (fabs(jets[i].eta()) < 5.0 && jets[i].pt() > 5.0) {
      n_good_jets++;
    }
  }
  results[14] = n_good_jets;

  // Z-Jet1 dEta and dR, only when |jeteta1|<5 and jetpt1>5
  results[15] = (jet_counter > 0 && fabs(jet1_eta) < 5.0 && jet1_pt > 5.0)
                    ? z_j1_eta
                    : -50.;
  results[16] = (jet_counter > 0 && fabs(jet1_eta) < 5.0 && jet1_pt > 5.0)
                    ? z_j1_dR
                    : -50.;

  // Jet1-Jet2 dR only when both jets pass cuts
  results[17] = (jet_counter > 1 && fabs(jet1_eta) < 5.0 && jet1_pt > 5.0 &&
                 fabs(jet2_eta) < 5.0 && jet2_pt > 5.0)
                    ? j1_j2_dR
                    : -50.;

  // Jet1-Jet3 and Jet2-Jet3 dR only when all three jets pass cuts
  results[18] = (jet_counter > 2 && fabs(jet1_eta) < 5.0 && jet1_pt > 5.0 &&
                 fabs(jet3_eta) < 5.0 && jet3_pt > 5.0)
                    ? j1_j3_dR
                    : -50.;
  results[19] = (jet_counter > 2 && fabs(jet2_eta) < 5.0 && jet2_pt > 5.0 &&
                 fabs(jet3_eta) < 5.0 && jet3_pt > 5.0)
                    ? j2_j3_dR
                    : -50.;

  return;
}