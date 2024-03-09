#include "shower.h"

Shower::Shower() {}

/**
 * This Function is Different from S.H's Tutorial, we keep it the same as the
 * GPU version for a fair test
 */
void Shower::SelectWinner(Event& ev, std::mt19937& gen) {
  std::uniform_real_distribution<> dis(0.0, 1.0);

  // Default Values
  double win_tt = tC;  // Lowest possible value is Cutoff Scale (in base.cuh)
  int win_sf = 16;     // 16 = No Splitting (0 -> 15 are Splitting Functions)
  int win_ij = 0;
  int win_k = 0;
  double win_zp = 0.0;
  double win_m2 = 0.0;

  for (int ij = 2; ij < ev.GetSize(); ij++) {
    for (int k = 2; k < ev.GetSize(); k++) {
      // Sanity Check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // Need to check if ij and k are colour connected
      if (!ev.GetParton(ij).IsColorConnected(ev.GetParton(k))) {
        continue;
      }

      // Params Identical to all splitting functions
      double m2 = (ev.GetParton(ij).GetMom() + ev.GetParton(k).GetMom()).M2();
      if (m2 < 4.0 * tC) {
        continue;
      }

      double zp = 0.5 * (1.0 + sqrt(1.0 - 4.0 * tC / m2));

      /**
       * If ij = quark/antiq, only one of the first 10 kernels can be used
       * If ij = gluon, we need to get one of six kernels
       */
      int sf = 16;  // Default or Null Splitting Function
      double tt = 0.0;

      // Different from S.H.'s Tutorial here, see shower.cu for more info.

      // Quark/Anti-Quark
      if (ev.GetParton(ij).GetPid() != 21) {
        // Flavour:  -5 -4 -3 -2 -1  1  2  3  4  5
        // Location:  0  1  2  3  4  5  6  7  8  9

        // Get the right splitting function
        sf = ev.GetParton(ij).GetPid() < 0 ? ev.GetParton(ij).GetPid() + 5
                                           : ev.GetParton(ij).GetPid() + 4;

        // Get t_temp
        double integral = kCF * 2.0 * std::log((1.0 - (1.0 - zp)) / (1.0 - zp));
        double g = asmax / (2.0 * M_PI) * integral;
        tt = ev.GetShowerT() * pow(dis(gen), 1.0 / g);
      } else if (ev.GetParton(ij).GetPid() == 21) {
        /**
         * We will need to give a chance to all six possible splittings
         *
         * 10: g -> gg
         * 11 to 15: g -> qqbar
         */

        double gg_integral = kCA * log((1.0 - (1.0 - zp)) / (1.0 - zp));
        double qq_integral = kTR / 2.0 * (zp - (1.0 - zp));

        double gg_g = asmax / (2.0 * M_PI) * gg_integral;
        double qq_g = asmax / (2.0 * M_PI) * qq_integral;

        // Get tt for g -> gg
        double gg_tt = ev.GetShowerT() * pow(dis(gen), 1.0 / gg_g);

        // Get tt for g -> qqbar
        // Select the quark flavour with highest random number
        double qq_tt = 0.0;
        int qnum = 0;
        for (int j = 1; j <= 5; j++) {
          // Get tt for g -> qqbar for the selected quark flavour
          double temp_tt = ev.GetShowerT() * pow(dis(gen), 1.0 / qq_g);

          if (temp_tt > qq_tt) {
            qq_tt = temp_tt;
            qnum = j;
          }
        }

        // Compare tt for g -> gg and g -> qqbar
        if (gg_tt > qq_tt) {
          sf = 10;
          tt = gg_tt;
        } else {
          sf = 10 + qnum;
          tt = qq_tt;
        }
      }

      // Check if tt is greater than the current winner
      if (sf < 16 && tt > win_tt) {
        win_tt = tt;
        win_sf = sf;
        win_ij = ij;
        win_k = k;
        win_zp = zp;
        win_m2 = m2;
      }
    }
  }

  // Store the results
  ev.SetShowerT(win_tt);
  ev.SetWinSF(win_sf);
  ev.SetWinDipole(0, win_ij);
  ev.SetWinDipole(1, win_k);
  ev.SetWinParam(0, win_zp);
  ev.SetWinParam(1, win_m2);
}

void Shower::MakeKinematics(Vec4* kinematics, const double z, const double y,
                            const double phi, const Vec4 pijt, const Vec4 pkt) {
  Vec4 Q = pijt + pkt;

  // Generating the Momentum (0, kt1, kt2, 0)
  double rkt = sqrt(Q.M2() * y * z * (1. - z));

  Vec4 kt1 = pijt.Cross(pkt);
  if (kt1.P() < 1.e-6) {
    Vec4 xaxis(0., 1., 0., 0.);
    kt1 = pijt.Cross(xaxis);
  }
  kt1 = kt1 * (rkt * cos(phi) / kt1.P());

  Vec4 kt2cms = Q.Boost(pijt);
  kt2cms = kt2cms.Cross(kt1);
  kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.P());
  Vec4 kt2 = Q.BoostBack(kt2cms);

  // Conversion to {i, j, k} basis
  Vec4 pi = pijt * z + pkt * ((1. - z) * y) + kt1 + kt2;
  Vec4 pj = pijt * (1. - z) + pkt * (z * y) - kt1 - kt2;
  Vec4 pk = pkt * (1. - y);

  // No need to do *kinematics[0], for arrays the elements are already pointers
  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
}

void Shower::MakeColours(Event& ev, int* coli, int* colj, const int flavs[3],
                         const int colij[2], const int colk[2], const int r) {
  // Increase variable ev.GetShowerC() by 1
  ev.IncrementShowerC();

  if (flavs[0] != 21) {
    if (flavs[0] > 0) {
      coli[0] = ev.GetShowerC();
      coli[1] = 0;
      colj[0] = colij[0];
      colj[1] = ev.GetShowerC();
    } else {
      coli[0] = 0;
      coli[1] = ev.GetShowerC();
      colj[0] = ev.GetShowerC();
      colj[1] = colij[1];
    }
  } else {
    if (flavs[1] == 21) {
      if (colij[0] == colk[1]) {
        if (colij[1] == colk[0] && r > 0.5) {
          coli[0] = colij[0];
          coli[1] = ev.GetShowerC();
          colj[0] = ev.GetShowerC();
          colj[1] = colij[1];
        } else {
          coli[0] = ev.GetShowerC();
          coli[1] = colij[1];
          colj[0] = colij[0];
          colj[1] = ev.GetShowerC();
        }
      } else {
        coli[0] = colij[0];
        coli[1] = ev.GetShowerC();
        colj[0] = ev.GetShowerC();
        colj[1] = colij[1];
      }
    } else {
      if (flavs[1] > 0) {
        coli[0] = colij[0];
        coli[1] = 0;
        colj[0] = 0;
        colj[1] = colij[1];
      } else {
        coli[0] = 0;
        coli[1] = colij[1];
        colj[0] = colij[0];
        colj[1] = 0;
      }
    }
  }
}

/**
 * In the GPU version, this would be split into multiple CUDA Kernels
 */
void Shower::GenerateSplitting(Event& ev, std::mt19937& gen) {
  std::uniform_real_distribution<> dis(0.0, 1.0);

  while (ev.GetShowerT() > tC) {
    SelectWinner(ev, gen);

    if (ev.GetShowerT() > tC) {
      double z = 0.0;
      double zp = ev.GetWinParam(0);
      double zm = 1.0 - zp;

      int kernel = ev.GetWinSF();
      int ij_flav = ev.GetParton(ev.GetWinDipole(0)).GetPid();

      double rand = dis(gen);
      if (kernel < 10 && ij_flav != 21) {
        z = 1.0 + (zp - 1.0) * pow((1.0 - zm) / (1.0 - zp), rand);
      } else if (kernel < 11) {
        z = 1.0 + (zp - 1.0) * pow((1.0 - zm) / (1.0 - zp), rand);
      } else if (kernel < 16) {
        z = zm + (zp - zm) * rand;
      }

      double y = ev.GetShowerT() / ev.GetWinParam(1) / z / (1.0 - z);

      double f = 0.0;
      double g = 0.0;
      double value = 0.0;
      double estimate = 0.0;

      // CS Kernel: y can't be 1
      if (y < 1.0) {
        if (kernel < 10 && ij_flav != 21) {
          value = kCF * (2.0 / (1.0 - z * (1.0 - y)) - (1.0 + z));
          estimate = kCF * 2.0 / (1.0 - z);
        } else if (kernel < 11) {
          value =
              kCA / 2.0 * (2.0 / (1.0 - z * (1.0 - y)) - 2.0 + z * (1.0 - z));
          estimate = kCA / (1.0 - z);
        } else if (kernel < 16) {
          value = kTR / 2.0 * (1.0 - 2.0 * z * (1.0 - z));
          estimate = kTR / 2.0;
        }

        f = (1.0 - y) * as(ev.GetShowerT()) * value;
        g = asmax * estimate;

        if (dis(gen) < f / g) {
          ev.SetShowerZ(z);
          ev.SetShowerY(y);

          double phi = 2.0 * M_PI * dis(gen);

          int win_ij = ev.GetWinDipole(0);
          int win_k = ev.GetWinDipole(1);

          Vec4 moms[3] = {Vec4(), Vec4(), Vec4()};
          MakeKinematics(moms, z, y, phi, ev.GetParton(win_ij).GetMom(),
                         ev.GetParton(win_k).GetMom());

          // Get Flavs from Kernel Number
          int kernel = ev.GetWinSF();

          // flavour: -5 -4 -3 -2 -1  1  2  3  4  5
          // index:    0  1  2  3  4  5  6  7  8  9
          int flavs[3];
          if (kernel < 10 && ev.GetParton(ev.GetWinDipole(0)).GetPid() != 21) {
            if (kernel < 5) {
              flavs[0] = kernel - 5;
              flavs[1] = kernel - 5;
              flavs[2] = 21;
            } else if (kernel < 10) {
              flavs[0] = kernel - 4;
              flavs[1] = kernel - 4;
              flavs[2] = 21;
            }
          } else if (kernel < 11) {
            flavs[0] = 21;
            flavs[1] = 21;
            flavs[2] = 21;
          } else if (kernel < 16) {
            flavs[0] = 21;
            flavs[1] = kernel - 10;
            flavs[2] = -1 * (kernel - 10);
          }
          int colij[2] = {ev.GetParton(win_ij).GetCol(),
                          ev.GetParton(win_ij).GetAntiCol()};

          int colk[2] = {ev.GetParton(win_k).GetCol(),
                         ev.GetParton(win_k).GetAntiCol()};

          int coli[2] = {0, 0};
          int colj[2] = {0, 0};
          MakeColours(ev, coli, colj, flavs, colij, colk, dis(gen));

          // Modify Splitter
          ev.SetPartonPid(win_ij, flavs[1]);
          ev.SetPartonMom(win_ij, moms[0]);
          ev.SetPartonCol(win_ij, coli[0]);
          ev.SetPartonAntiCol(win_ij, coli[1]);

          // Modify Recoiled Spectator
          ev.SetPartonMom(win_k, moms[2]);

          // Add Emitted Parton
          Parton em = Parton(flavs[2], moms[1], colj[0], colj[1]);
          ev.SetParton(ev.GetSize(), em);

          // Increment Emissions (IMPORTANT)
          ev.IncrementEmissions();

          return;
        }
      }
    }
  }
}

void Shower::Run(Event& ev) {

  /**
   * Thread Local
   * ------------
   *
   * We observed significant slowdown due to the rng - this is because the
   * code was re initialising the RNG. Using thread_local means that it is
   * initialised once and then reused, giving a massive speed-up!
   */
  thread_local std::random_device rd;
  thread_local std::mt19937 gen(rd());

  // Set the starting shower scale
  double t_max = (ev.GetParton(0).GetMom() + ev.GetParton(1).GetMom()).M2();
  ev.SetShowerT(t_max);

  // Set the initial number of emissions
  ev.SetEmissions(0);

  // Set the Colour Counter to 1 (q and qbar)
  ev.SetShowerC(1);

  while (ev.GetShowerT() > tC) {
    GenerateSplitting(ev, gen);
  }
}
