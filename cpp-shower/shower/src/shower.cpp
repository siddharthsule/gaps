#include "shower.h"

Shower::Shower() {}

/**
 * This Function is Different from S.H's Tutorial, we keep it the same as the
 * GPU version for a fair test
 */
void Shower::SelectWinner(Event& ev, std::mt19937& gen) {
  std::uniform_real_distribution<> dis(0., 1.);

  // Default Values
  double win_tt = tC;  // Lowest possible value is Cutoff Scale (in base.cuh)
  int win_sf = 0;      // 0 = No Splitting
  int win_ij = 0;
  int win_k = 0;
  double win_zp = 0.;
  double win_m2 = 0.;

  for (int ij = 0; ij < ev.GetSize(); ij++) {
    for (int k = 0; k < ev.GetSize(); k++) {
      // Sanity Check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // Skip non-partons
      if (!ev.GetParton(ij).IsParton() || !ev.GetParton(k).IsParton()) {
        continue;
      }

      // Need to check if ij and k are colour connected
      if (!ev.GetParton(ij).IsColorConnected(ev.GetParton(k))) {
        continue;
      }

      // Params Identical to all splitting functions
      double m2 = (ev.GetParton(ij).GetMom() + ev.GetParton(k).GetMom()).M2();
      if (m2 < 4. * tC) {
        continue;
      }

      double zp = 0.5 * (1. + sqrt(1. - 4. * tC / m2));

      // Codes instead of Object Oriented Approach!
      for (int sf : sfCodes) {
        // Check if the Splitting Function is valid for the current partons
        if (!validateSplitting(ev.GetParton(ij).GetPid(), sf)) {
          continue;
        }

        // Calculate the Evolution Variable
        double g = asmax / (2. * M_PI) * sfIntegral(1 - zp, zp, sf);
        double tt = ev.GetShowerT() * pow(dis(gen), 1. / g);

        // Check if tt is greater than the current winner
        if (tt > win_tt) {
          win_tt = tt;
          win_sf = sf;
          win_ij = ij;
          win_k = k;
          win_zp = zp;
          win_m2 = m2;
        }
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

/**
 * In the GPU version, this would be split into multiple CUDA Kernels
 */
void Shower::GenerateSplitting(Event& ev, std::mt19937& gen) {
  std::uniform_real_distribution<> dis(0., 1.);

  while (ev.GetShowerT() > tC) {
    SelectWinner(ev, gen);

    if (ev.GetShowerT() > tC) {
      // Get the Splitting Function
      int sf = ev.GetWinSF();

      double rand = dis(gen);

      // Generate z
      double zp = ev.GetWinParam(0);
      double z = sfGenerateZ(1 - zp, zp, rand, sf);

      double y = ev.GetShowerT() / ev.GetWinParam(1) / z / (1. - z);

      double f = 0.;
      double g = 0.;
      double value = 0.;
      double estimate = 0.;

      // CS Kernel: y can't be 1
      if (y < 1.) {
        value = sfValue(z, y, sf);
        estimate = sfEstimate(z, sf);

        f = (1. - y) * as(ev.GetShowerT()) * value;
        g = asmax * estimate;

        if (dis(gen) < f / g) {
          ev.SetShowerZ(z);
          ev.SetShowerY(y);

          double phi = 2. * M_PI * dis(gen);

          int win_ij = ev.GetWinDipole(0);
          int win_k = ev.GetWinDipole(1);

          Vec4 moms[3] = {Vec4(), Vec4(), Vec4()};
          MakeKinematics(moms, z, y, phi, ev.GetParton(win_ij).GetMom(),
                         ev.GetParton(win_k).GetMom());

          int flavs[3];
          sfToFlavs(sf, flavs);

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

void Shower::Run(Event& ev, int seed) {
  /**
   * Thread Local
   * ------------
   *
   * We observed significant slowdown due to the rng - this is because the
   * code was re initialising the RNG. Using thread_local means that it is
   * initialised once and then reused, giving a massive speed-up!
   */
  thread_local std::random_device rd;
  thread_local std::mt19937 gen(seed = -1 ? rd() : seed);

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
