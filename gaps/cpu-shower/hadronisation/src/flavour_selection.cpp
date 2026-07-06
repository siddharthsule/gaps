#include "hadronisation.h"

int hadronisation::select_qq_flavour(double rho, bool include_diquarks,
                                     double rho_diq) const {
  /**
   * @brief Select a flavour from the quark and diquark pools based on the
   * probability weights defined in the hadronisation class.
   *
   * we have pwt[0]=d, pwt[1]=u, pwt[2]=s and pwt[3]=diq. Our goal is to
   * select the flavour after normalising the weights
   *
   * pwt_total = pwt[0] + pwt[1] + pwt[2] + (include_diquarks ? pwt[3] : 0)
   *
   * rho is a random number between 0 and 1. We will use it to select the
   * flavour based on the cumulative distribution of the weights.
   *
   * @param include_diquarks whether to include diquarks in the selection pool
   * @param rho a random number between 0 and 1 for selection
   * @return the selected flavour (1=d, 2=u, 3=s, or diquark pid)
   */

  double pwt_tot = pwt[0] + pwt[1] + pwt[2] + (include_diquarks ? pwt[3] : 0);
  double r = rho * pwt_tot;

  // d
  if (r < pwt[0]) {
    return 1;
  }

  // u
  if (r < pwt[0] + pwt[1]) {
    return 2;
  }

  // s
  if (r < pwt[0] + pwt[1] + pwt[2]) {
    return 3;
  }

  // diquarks
  if (include_diquarks) {
    // Map rho_diq to the diquark pool (all light diquarks, once each)
    const int diquark_pool[9] = {
        1103,         // (dd)_1  vector
        2101, 2103,   // (ud)_0  scalar, (ud)_1  vector
        2203,         // (uu)_1  vector
        3101, 3103,   // (sd)_0  scalar, (sd)_1  vector
        3201, 3203,   // (su)_0  scalar, (su)_1  vector
        3303          // (ss)_1  vector
    };

    int r2 = static_cast<int>(rho_diq * 9.0);
    if (r2 > 8) r2 = 8;
    return diquark_pool[r2];
  }

  // Fallback (should not reach here)
  return 0;
}