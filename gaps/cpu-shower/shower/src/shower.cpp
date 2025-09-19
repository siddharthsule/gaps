#include "shower.h"

// -----------------------------------------------------------------------------

void shower::select_winner(event& ev, double* winner) const {
  /**
   * @brief Select the winner splitting in the event
   *
   * This function generates the highest transverse momentum splitting for every
   * possible dipoles in the event. It then chooses the winner emission from the
   * generated splittings, by picking the one with the highest transverse
   * momentum. This winner emission is then used in the veto step
   *
   * @param ev The event to select the winner from
   * @param winner The array to store the winner emission data
   */

  // default values
  double win_tt = t_c;  // lowest possible value is cutoff scale (in base.cuh)
  int win_sf = 0;       // 0 = no splitting
  int win_ij = 0;
  int win_k = 0;
  double win_sijk = 0.;
  double win_zm = 0.;
  double win_zp = 0.;

  for (int ij = 0; ij < ev.get_size(); ij++) {
    for (int k = 0; k < ev.get_size(); k++) {
      // sanity check to ensure ij != k
      if (ij == k) {
        continue;
      }

      // Skip non-partons
      if (!ev.get_particle(ij).is_parton() || !ev.get_particle(k).is_parton()) {
        continue;
      }

      // need to check if ij and k are colour connected
      if (!ev.get_particle(ij).is_color_connected(ev.get_particle(k))) {
        continue;
      }

      // params identical to all splitting functions
      double sijk =
          (ev.get_particle(ij).get_mom() + ev.get_particle(k).get_mom()).m2();
      if (sijk < 4. * t_c) {
        continue;
      }

      // get the splitting functions for the current partons
      // sf_codes is an array of possible splitting functions
      int sf_codes[11];
      generate_possible_splittings(ev.get_particle(ij).get_pid(),
                                   ev.get_particle(k).get_pid(),
                                   ev.get_particle(ij).is_initial(),
                                   ev.get_particle(k).is_initial(), sf_codes);

      // codes instead of object oriented approach! See splittings.cpp
      for (int sf : sf_codes) {
        // When a null code is encountered, we have reached the end of the
        // possible splittings, we can break out of the loop
        if (sf == -1) {
          break;
        }

        // check if either parton has eta < eta_min (usually 1e-5)
        if (is_ii(sf) && (ev.get_particle(ij).get_eta() < 1e-5 ||
                          ev.get_particle(k).get_eta() < 1e-5)) {
          continue;
        }

        // phase space limits
        double zm, zp;
        get_boundaries(zm, zp, sijk, ev.get_particle(ij).get_eta(), sf);

        if (zm < 0. || zp > 1. || zm > zp) {
          continue;
        }

        // Calulate the integrated overestimate
        double pdf_max = get_pdf_max(sf, ev.get_particle(ij).get_eta());
        double c =
            as_max / (2. * M_PI) * sf_integral(zm, zp, sf) * j0_max * pdf_max;

        // calculate the evolution variable
        double tt;

        // Final State Radiation - As always
        if (is_fsr(sf)) {
          // t = T * random^(1/c)
          tt = ev.get_shower_t() * pow(ev.gen_random(), 1. / c);
        }

        // Initial State Radiation - Need to account for quark masses for PDF
        // Ranges (Avoid evolving to below quark mass)
        else {
          // Get the quark mass squared
          int fl = abs(ev.get_particle(ij).get_pid());
          double mq2 = fl == 5 ? mb * mb : (fl == 4 ? mc * mc : 0.);

          // Evolve from (t - m2) not t for this case
          tt = (ev.get_shower_t() - mq2) * pow(ev.gen_random(), 1. / c) + mq2;

          // Check if tt <= mq2: Need some numerical error protection
          if (tt - mq2 < 1e-9) {
            continue;
          }
        }

        // For FI, IF and II, skip events where q2 (=tt) is less than pdf limit
        if (!is_ff(sf) && tt < pdf_q_min * pdf_q_min) {
          continue;
        }

        // check if tt is greater than the current winner
        if (tt > win_tt) {
          win_tt = tt;
          win_sf = sf;
          win_ij = ij;
          win_k = k;
          win_sijk = sijk;
          win_zm = zm;
          win_zp = zp;
        }
      }
    }
  }

  // set the new shower t
  ev.set_shower_t(win_tt);

  // Also generate z, y and phi
  double z = sf_generate_z(win_zm, win_zp, ev.gen_random(), win_sf);
  double y = calculate_y(win_tt, z, win_sijk, win_sf);
  double phi = 2. * M_PI * ev.gen_random();

  // IF: z, y -> x, u (here, z, y)
  if (is_if(win_sf)) {
    double z0 = z;
    double ratio = win_tt / win_sijk;
    double frac2 =
        4 * ratio * z0 * (1. - z0) / ((1. - z0 + ratio) * (1. - z0 + ratio));
    z = (1. - z0 + ratio) / (2 * ratio) * (1. - sqrt(1. - frac2));
    y = (z * ratio) / (1. - z0);
  }

  // II: z, pt -> x, v (here, z, y)
  else if (is_ii(win_sf)) {
    double z0 = z;
    double ratio = win_tt / win_sijk;
    z = z0 * (1. - z0) / (1. - z0 + ratio);
    y = (z * ratio) / (1. - z0);
  }

  // Set the winner variables (sf, ij, k, sijk, z, y, phi)
  winner[0] = static_cast<double>(win_sf);
  winner[1] = static_cast<double>(win_ij);
  winner[2] = static_cast<double>(win_k);
  winner[3] = win_sijk;
  winner[4] = z;
  winner[5] = y;
  winner[6] = phi;

  return;
}

// -----------------------------------------------------------------------------

/**
 * in the gpu version, this would be split into multiple cuda kernels
 */

void shower::generate_splitting(event& ev) {
  /**
   * @brief Generate a splitting for the current event
   *
   * This function generates a splitting for the current event, by selecting the
   * winner emission and performing the veto algorithm. If the veto is passed,
   * the event is modified to include the new parton and the shower scale is
   * updated.
   *
   * @param ev The event to generate the splitting for
   */

  while (ev.get_shower_t() > t_c) {
    /**
     * Shower Variables - useful to store as collective
     *
     * t, c and end_shower stored in event, becuase they
     * are unique to each event, and not throwaway values
     * like these.
     *
     * Winner variables: (sf, ij, k, sijk, z, y, phi)
     * Stored in ONE array, so we make it 6 x n_events
     * Stored all as doubles, so static_cast<int> for sf, ij, k
     */
    double winner[7] = {0., 0., 0., 0., 0., 0., 0.};

    // select the winner splitting function
    select_winner(ev, winner);

    // Get the shower evolution variable
    double t = ev.get_shower_t();

    // Get the winner variables (sf, ij, j, sijk, z, y, phi)
    int sf = static_cast<int>(winner[0]);
    int ij = static_cast<int>(winner[1]);
    int k = static_cast<int>(winner[2]);
    double sijk = winner[3];
    double z = winner[4];
    double y = winner[5];
    double phi = winner[6];

    // calculate z (as in the sampled variable)
    double z0;
    if (is_isr(sf)) {
      z0 = 1. - (z * (t / sijk) / y);
    } else {
      // FF and FI: can use z directly
      z0 = z;
    }

    // Check for Cutoff
    if (!(t > t_c)) {
      // set the end shower flag
      ev.shower_has_ended(true);
      break;
    }

    // Check for Phase Space
    if (!(check_phase_space(z, y, sf))) {
      continue;
    }

    // Additional Check for FI and IF/II, in case quark > proton
    vec4 pijt = ev.get_particle(ij).get_mom();
    vec4 pkt = ev.get_particle(k).get_mom();
    if ((is_isr(sf)) && (pijt[0] / z > e_proton)) {
      continue;
    }
    if ((is_fi(sf)) && (pkt[0] / y > e_proton)) {
      continue;
    }

    // Get PDF Ratio and PDF Max for FI and IF/II
    double pdf_ratio(1.), pdf_max(1.);

    if (!is_ff(sf)) {
      // Check Momentum Fraction
      if (!(check_mom_frac(sf, ev.get_particle(ij).get_pid(),
                           ev.get_particle(k).get_pid(),
                           ev.get_particle(ij).get_eta(),
                           ev.get_particle(k).get_eta(), z))) {
        continue;
      }

      // Get PDF Ratio
      pdf_ratio = get_pdf_ratio(
          sf, ev.get_particle(ij).get_pid(), ev.get_particle(k).get_pid(),
          ev.get_particle(ij).get_eta(), ev.get_particle(k).get_eta(), z, t);

      // cancel emission if pdf_ratio is nan, inf, or less than 0
      if (isnan(pdf_ratio) || isinf(pdf_ratio) || pdf_ratio <= 0.) {
        continue;
      }

      // If PDF Ratio is too large, veto
      if (pdf_ratio > 1e6) {
        continue;
      }

      /**
       * Why the factor of z?
       * --------------------
       * From LHAPDF, we get xf(x, q2). This means our ratio will be equal to
       * (x/z)f(x/z, q2) / xf(x, q2) = 1/z * f(x/z, q2) / f(x, q2)
       *
       * So, we can either adjust the jacobian by multiplying by z or adjust the
       * pdf_ratio by dividing by z. We choose the latter.
       */
      pdf_ratio *= z;

      // Mutliply by (t - m2) / t for ISR to account for quark masses
      int fl = abs(ev.get_particle(ij).get_pid());
      pdf_ratio *= (t - (fl == 5 ? mb * mb : (fl == 4 ? mc * mc : 0.))) / t;

      // -----------------------------------------------------------------------
      /*
      // Printing for PDF Ratio Analysis
      int pid_final;
      if (get_splitting_type(sf) == 0) {
        pid_final = ev.get_particle(ij).get_pid();
      } else if (get_splitting_type(sf) == 1) {
        pid_final = 21;
      } else if (get_splitting_type(sf) == 2) {
        pid_final = ev.get_particle(ij).get_pid();
      } else if (get_splitting_type(sf) == 3) {
        pid_final = -1 * get_splitting_flavour(sf);
      } else if (get_splitting_type(sf) == 4) {
        pid_final = get_splitting_flavour(sf);
      }

      double eta_out;
      if (is_fi(sf) == 1) {
        eta_out = ev.get_particle(k).get_eta();
      } else {
        eta_out = ev.get_particle(ij).get_eta();
      }

      std::cout << ev.get_particle(ij).get_pid() << ", " << pid_final << ", "
                << eta_out << ", " << pdf_ratio << std::endl;
      */
      // -----------------------------------------------------------------------

      // Get PDF Max
      pdf_max = get_pdf_max(sf, ev.get_particle(ij).get_eta());
    }

    // Jacobian
    double jacobian = get_jacobian(z, y, z0, sf) * pdf_ratio;
    double jmaxtot = j0_max * pdf_max;

    // Splitting Function Value and Estimate
    double value = sf_value(z, y, sf);
    double estimate = sf_estimate(z0, sf);

    // veto algorithm
    double f = as(t) * value * jacobian;
    double g = as_max * estimate * jmaxtot;

    // Check for Negative f
    if (f < 0.) {
      continue;
    }

    // Check for f > g
    if (f > g) {
      continue;
    }

    // Veto Algorithm
    if (ev.gen_random() < f / g) {
      // Success!

      // get the flavours
      int flavs[3];
      sf_to_flavs(sf, flavs);

      // pi, pj, pk, pijt, pkt and kt
      vec4 moms[6] = {vec4(), vec4(), vec4(), vec4(), vec4(), vec4()};
      make_kinematics(moms, z, y, phi, ev.get_particle(ij).get_mom(),
                      ev.get_particle(k).get_mom(), sf);

      // calculate the colours
      int colij[2] = {ev.get_particle(ij).get_col(),
                      ev.get_particle(ij).get_acol()};
      int colk[2] = {ev.get_particle(k).get_col(),
                     ev.get_particle(k).get_acol()};
      int coli[2] = {0, 0};
      int colj[2] = {0, 0};
      make_colours(ev.get_shower_c(), sf, flavs, colij, colk, coli, colj,
                   ev.gen_random());

      // modify splitter
      ev.set_particle_pid(ij, flavs[1]);
      ev.set_particle_mom(ij, moms[0]);
      ev.set_particle_col(ij, coli[0]);
      ev.set_particle_acol(ij, coli[1]);
      if (is_isr(sf)) {
        ev.set_particle_eta(ij, ev.get_particle(ij).get_eta() / z);
      }

      // modify recoiled spectator
      ev.set_particle_mom(k, moms[2]);
      if (is_fi(sf)) {
        ev.set_particle_eta(k, ev.get_particle(k).get_eta() / y);
      }

      // add emitted parton
      particle em = particle(flavs[2], moms[1], colj[0], colj[1]);
      ev.set_particle(ev.get_size(), em);

      // increment emissions (important) !!!!!!
      ev.increment_emissions();

      // II Only - Lorentz Boost the new final state
      if (is_ii(sf)) {
        ii_boost_after_emission(ev, moms);
      }

      return;
    }
  }
}

void shower::run(event& ev, bool nlo_matching) {
  /**
   * @brief Run the shower for the current event
   *
   * This function runs the shower for the current event, by generating
   * splittings until the shower scale is below the cutoff scale. The shower
   * scale is set to the smallest transverse momentum splitting in the event,
   * in order to prevent phase space overlap.
   *
   * @param ev The event to run the shower for
   */

  // NLO Matching does the first emission and sets the shower scale
  // to the first emission pT. If NLO Matching is off, we find the
  // smallest pT in the event and set the shower scale to that.
  if (!nlo_matching) {
    // set the starting shower scale
    double t_start = 10000000.;
    for (int i = 0; i < ev.get_size(); i++) {
      for (int j = 0; j < ev.get_size(); j++) {
        if (i == j) {
          continue;
        }

        if (!ev.get_particle(i).is_parton() ||
            !ev.get_particle(j).is_parton()) {
          continue;
        }

        double t =
            (ev.get_particle(i).get_mom() + ev.get_particle(j).get_mom()).m2();

        if (t < t_start) {
          t_start = t;
        }
      }
    }

    ev.set_shower_t(t_start);
    ev.set_shower_c(1);
  }

  // run the shower
  while (ev.get_shower_t() > t_c) {
    // limit to max particles
    if (ev.get_size() == min(max_particles, ev.get_hard() + n_emissions_max)) {
      // Only print warning if too many emissions for the code,
      // not when the number of emissions is limited by the user
      if (max_particles < n_emissions_max) {
        std::cout << "Warning: Max Particles Reached" << std::endl;
      }
      break;
    }

    generate_splitting(ev);
  }
}