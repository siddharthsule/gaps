#include "shower.h"

shower::shower() {}

/**
 * this function is different from s.h's tutorial, we keep it the same as the
 * gpu version for a fair test
 */
void shower::select_winner(event& ev) {
  // Seed and Random
  unsigned long seed = ev.get_seed();
  double rand = ev.get_rand();

  // default values
  double win_tt = t_c;  // lowest possible value is cutoff scale (in base.cuh)
  int win_sf = 0;       // 0 = no splitting
  int win_ij = 0;
  int win_k = 0;
  double win_zp = 0.;
  double win_m2 = 0.;

  // we start at 2 because elements 0 and 1 are electrons - to change with isr
  for (int ij = 2; ij < ev.get_size(); ij++) {
    for (int k = 2; k < ev.get_size(); k++) {
      // sanity check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // need to check if ij and k are colour connected
      if (!ev.get_parton(ij).is_color_connected(ev.get_parton(k))) {
        continue;
      }

      // params identical to all splitting functions
      double m2 =
          (ev.get_parton(ij).get_mom() + ev.get_parton(k).get_mom()).m2();
      if (m2 < 4. * t_c) {
        continue;
      }

      double zp = 0.5 * (1. + sqrt(1. - 4. * t_c / m2));

      // codes instead of object oriented approach!
      for (int sf : sf_codes) {
        // check if the splitting function is valid for the current partons
        if (!validate_splitting(ev.get_parton(ij).get_pid(), sf)) {
          continue;
        }

        // calculate the evolution variable
        double g = asmax / (2. * M_PI) * sf_integral(1 - zp, zp, sf);
        double tt = ev.get_shower_t() * pow(rand, 1. / g);
        seed = ev.get_seed();
        rand = ev.get_rand();
        update_rng(seed, rand);
        ev.set_seed(seed);
        ev.set_rand(rand);
        printf("Seed: %lu, Rand: %f, Select Winner \n", seed, rand);

        // check if tt is greater than the current winner
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

  // store the results
  ev.set_shower_t(win_tt);
  ev.set_win_sf(win_sf);
  ev.set_win_dipole(0, win_ij);
  ev.set_win_dipole(1, win_k);
  ev.set_win_param(0, win_zp);
  ev.set_win_param(1, win_m2);

  // set the seed and random number
  ev.set_seed(seed);
  ev.set_rand(rand);
}

/**
 * in the gpu version, this would be split into multiple cuda kernels
 */
void shower::generate_splitting(event& ev) {
  // Seed and Random
  unsigned long seed = ev.get_seed();
  double rand = ev.get_rand();

  while (ev.get_shower_t() > t_c) {
    select_winner(ev);

    if (ev.get_shower_t() > t_c) {
      // get the splitting function
      int sf = ev.get_win_sf();

      // generate z
      double zp = ev.get_win_param(0);
      double z = sf_generate_z(1 - zp, zp, rand, sf);
      seed = ev.get_seed();
      rand = ev.get_rand();
      update_rng(seed, rand);
      ev.set_seed(seed);
      ev.set_rand(rand);
      printf("Seed: %lu, Rand: %f, Z Gen\n", seed, rand);

      double y = ev.get_shower_t() / ev.get_win_param(1) / z / (1. - z);

      double f = 0.;
      double g = 0.;
      double value = 0.;
      double estimate = 0.;

      // cs kernel: y can't be 1
      if (y < 1.) {
        value = sf_value(z, y, sf);
        estimate = sf_estimate(z, sf);

        f = (1. - y) * as(ev.get_shower_t()) * value;
        g = asmax * estimate;

        // Just to Avoid Confusion...
        double r = rand;
        seed = ev.get_seed();
        rand = ev.get_rand();
        update_rng(seed, rand);
        ev.set_seed(seed);
        ev.set_rand(rand);
        printf("Seed: %lu, Rand: %f, Veto Step\n", seed, rand);

        if (r < f / g) {
          ev.set_shower_z(z);
          ev.set_shower_y(y);

          double phi = 2. * M_PI * rand;
          seed = ev.get_seed();
          rand = ev.get_rand();
          update_rng(seed, rand);
          ev.set_seed(seed);
          ev.set_rand(rand);
          printf("Seed: %lu, Rand: %f, Phi Shower Gen\n", seed, rand);

          int win_ij = ev.get_win_dipole(0);
          int win_k = ev.get_win_dipole(1);

          vec4 moms[3] = {vec4(), vec4(), vec4()};
          make_kinematics(moms, z, y, phi, ev.get_parton(win_ij).get_mom(),
                          ev.get_parton(win_k).get_mom());

          int flavs[3];
          sf_to_flavs(sf, flavs);

          int colij[2] = {ev.get_parton(win_ij).get_col(),
                          ev.get_parton(win_ij).get_anti_col()};

          int colk[2] = {ev.get_parton(win_k).get_col(),
                         ev.get_parton(win_k).get_anti_col()};

          int coli[2] = {0, 0};
          int colj[2] = {0, 0};
          make_colours(ev, coli, colj, flavs, colij, colk, rand);
          seed = ev.get_seed();
          rand = ev.get_rand();
          update_rng(seed, rand);
          ev.set_seed(seed);
          ev.set_rand(rand);
          printf("Seed: %lu, Rand: %f, Colour Gen\n", seed, rand);

          // modify splitter
          ev.set_parton_pid(win_ij, flavs[1]);
          ev.set_parton_mom(win_ij, moms[0]);
          ev.set_parton_col(win_ij, coli[0]);
          ev.set_parton_anti_col(win_ij, coli[1]);

          // modify recoiled spectator
          ev.set_parton_mom(win_k, moms[2]);

          // add emitted parton
          parton em = parton(flavs[2], moms[1], colj[0], colj[1]);
          ev.set_parton(ev.get_size(), em);

          // increment emissions (important)
          ev.increment_emissions();

          // set the seed and random number - here if it returns early
          ev.set_seed(seed);
          ev.set_rand(rand);

          return;
        }
      }
    }
  }

  // set the seed and random number - if no emissions
  ev.set_seed(seed);
  ev.set_rand(rand);
}

void shower::run(event& ev) {
  // same seed option. turn off by commenting when not in use!
  // having an if statement if no seed is given would not be a fair comparison
  // to the gpu, so commented out is better for now. maybe in the future.
  // thread_local std::mt19937 gen(seed);

  // set the starting shower scale
  double t_max = (ev.get_parton(0).get_mom() + ev.get_parton(1).get_mom()).m2();
  ev.set_shower_t(t_max);

  // set the initial number of emissions
  ev.set_emissions(0);

  // set the colour counter to 1 (q and qbar)
  ev.set_shower_c(1);

  while (ev.get_shower_t() > t_c) {
    generate_splitting(ev);
  }
}
