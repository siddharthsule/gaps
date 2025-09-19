#include "shower.cuh"

__device__ void shower::get_boundaries(double& zm, double& zp, double sijk,
                                       double eta, int sf) const {
  /**
   * @brief Get the phase space boundaries for the given splitting
   *
   * @param zm: Lower boundary
   * @param zp: Upper boundary
   * @param sijk: Invariant mass of the emitter and spectator
   * @param eta: Momentum fraction of the emitter
   * @param sf: Splitting function code
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  // eta_max is the max x value of PDF set
  // double eta_max = pdf_x_max;

  double t_max, frac;

  switch (splitting_case) {
    // FF splittings: z and y
    case 0:
      zp = 0.5 * (1. + sqrt(1. - 4. * t_c / sijk));
      zm = 1. - zp;
      return;
      break;

    // FI splittings: z and x (here, y)
    case 1:
      zp = 0.5 * (1. + sqrt(1. - 4. * t_c / sijk));
      zm = 1. - zp;
      return;
      break;

    // IF splittings: x (here, z) and u (here, y)
    case 2:
      // These are x boundaries
      // zm = eta / eta_max;
      // zp = (sijk) / (sijk + 4 * t_c);

      // These are z boundaries
      frac = (4 * t_c * eta) / (sijk * (1. - eta));
      zm = 0.5 * (1. + eta - (1. - eta) * sqrt(1. - frac));
      zp = 0.5 * (1. + eta + (1. - eta) * sqrt(1. - frac));
      return;
      break;

    // II splittings: x (here, z) and v (here, y)
    case 3:
      // These are x boundaries
      // zm = eta / eta_max;
      // zp = (sijk) / (sijk + 4 * t_c);

      // These are z boundaries
      t_max = (sijk * (1. - eta) * (1. - eta)) / (4 * eta);
      zm = 0.5 * (1. + eta - (1. - eta) * sqrt(1. - (t_c / t_max)));
      zp = 0.5 * (1. + eta + (1. - eta) * sqrt(1. - (t_c / t_max)));
      return;
      break;
  }
  zm = 0.;
  zp = 1.;
  return;
}

__device__ bool shower::check_phase_space(double z, double y, int sf) const {
  /**
   * @brief Check if the phase space boundaries are satisfied
   *
   * @param z: Splitting variable
   * @param y: Momentum fraction of the spectator
   * @param sf: Splitting function code
   * @return true: If the phase space boundaries are satisfied
   */

  // additional sanity check: check if y is not nan
  if (isnan(y)) {
    return false;
  }

  // luckily, all four splitting types have the same phase space boundaries
  if (z < 0. || z > 1. || y < 0. || y > 1.) {
    return false;
  }

  if (is_ii(sf)) {
    if (y > 1. - z) {
      return false;
    }
  }

  return true;
}