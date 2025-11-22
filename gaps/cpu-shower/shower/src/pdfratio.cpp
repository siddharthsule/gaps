#include "shower.h"

// -----------------------------------------------------------------------------
// Calculate the pdf ratio

double shower::pdf_ratio(int pid_before, int pid_after, double eta, double z,
                         double t) const {
  /**
   * @brief Calculate the ratio of the pdfs for pre and post emission particles
   *
   * @param pid_before Particle ID of the post-emission particle
   * @param pid_after Particle ID of the pre-emission particle
   * @param eta Momentum Fraction of the pre-emission particle
   * @param z Splitting Variable
   * @param t Evolution Variable
   * @return double: PDF Ratio
   */

  double xf_a = pdf->xfxQ2(pid_after, eta / z, t);  // After Emission
  double xf_b = pdf->xfxQ2(pid_before, eta, t);     // Before Emission

  // Ensure PDFs are valid
  if (isnan(xf_a) || isnan(xf_b) || isinf(xf_a) || isinf(xf_b) || xf_a <= 0. ||
      xf_b <= 0.) {
    return -1.;
  }

  // printf("PDFs: %f %f %f %f\n", xf_a, xf_b, z, t);
  return xf_a / xf_b;
}

// -----------------------------------------------------------------------------
// Wrapper to get the right pdf ratio

double shower::get_pdf_ratio(int sf, int ij_pid, int k_pid, double ij_eta,
                             double k_eta, double z, double t) const {
  /**
   * @brief Get the appropriate pdf ratio for the given splitting
   *
   * @param sf Splitting Function Code
   * @param ij_pid Particle ID of the emitter
   * @param k_pid Particle ID of the spectator
   * @param ij_eta Momentum Fraction of the emitter
   * @param k_eta Momentum Fraction of the spectator
   * @param z Splitting Variable
   * @param t Evolution Variable
   * @return double: PDF Ratio
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  // 0 = q->qg, 1 = q->gq, 2 = g->gg, 3 = g->qq
  int splitting_type = get_splitting_type(sf);

  // 0 = particle, 1 = antiparticle
  int is_em_antip = is_emitter_antiparticle(sf);

  // 1 = d, 2 = u, 3 = s, 4 = c, 5 = b; 0 = gluon
  int splitting_flavour = get_splitting_flavour(sf);

  switch (splitting_case) {
    // FI splittings: z and x (here, y)
    case 1:
      return pdf_ratio(k_pid, k_pid, k_eta, z, t);
      break;

    // IF splittings: x (here, z) and u (here, y)
    // II splittings: x (here, z) and v (here, y)
    case 2:
    case 3:
      // q -> qg
      if (splitting_type == 0) {
        return pdf_ratio(ij_pid, ij_pid, ij_eta, z, t);
      }

      // q -> gq (back wards is g -> q qbar / qbar q)
      else if (splitting_type == 1) {
        return pdf_ratio(ij_pid, 21, ij_eta, z, t);
      }

      // g -> gg
      else if (splitting_type == 2) {
        return pdf_ratio(ij_pid, ij_pid, ij_eta, z, t);
      }

      // g -> q qbar (backwards is qbar to g qbar)
      else if (splitting_type == 3) {
        return pdf_ratio(21, -1 * splitting_flavour, ij_eta, z, t);
      }

      else if (splitting_type == 4) {
        // g -> qbar q (backwards is q -> g q)
        return pdf_ratio(21, splitting_flavour, ij_eta, z, t);
      }

      break;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Check if the momentum fraction is valid before calculating

bool shower::check_mom_frac(int sf, int ij_pid, int k_pid, double ij_eta,
                            double k_eta, double z) const {
  /**
   * @brief Check if the momentum fraction is valid before calculating the pdf
   * ratio
   *
   * @param sf Splitting Function Code
   * @param ij_pid Particle ID of the emitter
   * @param k_pid Particle ID of the spectator
   * @param ij_eta Momentum Fraction of the emitter
   * @param k_eta Momentum Fraction of the spectator
   * @param z Splitting Variable
   * @return true: If the momentum fraction is valid
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  switch (splitting_case) {
    // FI splittings: z and x (here, y)
    case 1:
      return (k_eta < z) && (k_eta > pdf_x_min) && (k_eta < pdf_x_max) &&
             (k_eta / z > pdf_x_min) && (k_eta / z < pdf_x_max);
      break;

    // IF splittings: x (here, z) and u (here, y)
    case 2:
    case 3:
      return (ij_eta < z) && (ij_eta > pdf_x_min) && (ij_eta < pdf_x_max) &&
             (ij_eta / z > pdf_x_min) && (ij_eta / z < pdf_x_max);
      break;
  }
  return false;
}

// -----------------------------------------------------------------------------
// Get the appropriate pdf_max value

double shower::get_pdf_max(int sf, double ij_eta) const {
  /**
   * @brief Get the maximum pdf ratio for the given splitting
   *
   * @param sf Splitting Function Code
   * @param ij_eta Momentum Fraction of the emitter
   * @return double: PDF Max
   */

  // No pdf_max for FF splittings
  if (get_splitting_case(sf) == 0) {
    return 1.;
  }

  // ---------------------------------------------------------------------------
  // Extra Large PDF Maxes for PDF Ratio Analysis
  /*
  // IF/II q -> gq type splittings
  if ((get_splitting_case(sf) > 1) && (get_splitting_type(sf) == 1)) {
    switch (get_splitting_flavour(sf)) {
      case 5:
      case 4:
      case 3:
        return 10000.;
      default:
        return 1000.;
    }
  }

  else {
    return 10.;
  }
  */
  // ---------------------------------------------------------------------------

  /**
   * You can attain the following values for the pdf_max by running the
   * commented code above and in shower.cpp. For more information, refer
   * to the paper or the docs.
   */

  // IF/II q -> gq type splittings
  if ((get_splitting_case(sf) > 1) && (get_splitting_type(sf) == 1)) {
    // qbar -> g qbar
    if (is_emitter_antiparticle(sf)) {
      switch (get_splitting_flavour(sf)) {
        case 5:
          return max(150., 700. * ij_eta);
        case 4:
          return max(100., 450. * ij_eta);
        case 3:
          return max(50., 600. * ij_eta);
        default:
          return 150.;
      }
    }

    // q -> g q
    else {
      switch (get_splitting_flavour(sf)) {
        case 5:
          return max(150., 650. * ij_eta);
        case 4:
          return max(100., 450. * ij_eta);
        case 3:
          return max(50., 600. * ij_eta);
        default:
          return 150.;
      }
    }
  }

  // IF/II g -> d dbar
  else if ((sf == 2301) || (sf == 3301)) {
    return 5. * sqrt(ij_eta);
  }

  // IF/II g -> u ubar
  else if ((sf == 2302) || (sf == 3302)) {
    return 10. * sqrt(ij_eta);
  }

  // All other Splittings, y99 ~ y9999 ~ 1
  else {
    return 1.;
  }
}