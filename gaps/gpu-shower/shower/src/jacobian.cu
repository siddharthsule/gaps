#include "shower.cuh"

__device__ double shower::get_jacobian(double z, double y, int sf) const {
  /**
   * @brief Calculate the jacobian for the given splitting
   *
   * @param z: Splitting variable
   * @param y: Momentum fraction of the spectator
   * @param sf: Splitting function code
   * @return double: Jacobian
   */

  // 0 = FF, 1 = FI, 2 = IF, 3 = II
  int splitting_case = get_splitting_case(sf);

  switch (splitting_case) {
    // FF splittings: z and y
    case 0:
      return 1. - y;
      break;

    // FI splittings: z and x (here, y)
    case 1:
      return 1.;
      break;

    // IF splittings: x (here, z) and u (here, y)
    case 2:
      return 1. / (z + y - 2. * z * y);
      break;

    // II splittings: x (here, z) and v (here, y)
    case 3:
      return 1. / (z + y);
      break;
  }
  return 0.;
}