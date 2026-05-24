#ifndef qcd_h_
#define qcd_h_

#include "base.h"

// QCD constants
const double k_tr = 0.5;
const double k_ca = k_nc;
const double k_cf = (k_nc * k_nc - 1.) / (2. * k_nc);

class alpha_s {
  /**
   * @class alpha_s
   * @brief The strong coupling constant
   *
   * This class is used to calculate the strong coupling constant at a given
   * scale.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  int n_loops;
  bool use_cmw;
  double asmz, asmb, asmc;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  alpha_s(double asmz, int n_loops = 2, bool use_cmw = false)
      : n_loops(n_loops), use_cmw(use_cmw), asmz(asmz), asmb(0.), asmc(0.) {
    /**
     * @brief Construct the alpha_s object and pre-compute threshold values
     *
     * @param asmz alpha_s at the Z mass (MS-bar)
     * @param n_loops number of loops for the running coupling (0=fixed, 1=LO,
     * 2=NLO)
     * @param use_cmw apply the CMW scheme rescaling when evaluating alpha_s
     *
     * asmb and asmc are computed WITHOUT the CMW correction because CMW is
     * an additive shift applied at evaluation time, not a change to the
     * running itself. The threshold values must stay in the MS-bar scheme
     * so that crossing mb or mc gives a continuous coupling.
     */
    switch (n_loops) {
      case 0:
        asmb = asmz;
        asmc = asmz;
        break;
      case 1:
        asmb = as1(mb2);
        asmc = as1(mc2);
        break;
      case 2:
      default:
        asmb = as2(mb2);
        asmc = as2(mc2);
        break;
    }
  }

  // ---------------------------------------------------------------------------
  // member functions

  double beta0(int nf) const {
    /**
     * @brief calculate the beta function at 1 loop (leading order)
     *
     * @param nf the number of flavours
     * @return the beta function
     */

    return (11. / 6. * k_ca) - (2. / 3. * k_tr * nf);
  }

  double beta1(int nf) const {
    /**
     * @brief calculate the beta function at 2 loops (next-to-leading order)
     *
     * @param nf the number of flavours
     * @return the beta function
     */

    return (17. / 6. * k_ca * k_ca) - ((5. / 3. * k_ca + k_cf) * k_tr * nf);
  }

  double as1(double t) const {
    /**
     * @brief calculate the strong coupling constant at 1 loop (leading order)
     *
     * @param t the scale
     * @return the strong coupling constant
     */

    double tref, asref, b0;

    if (t >= mb2) {
      tref = mz2;
      asref = asmz;
      b0 = beta0(5) / (2. * M_PI);
    } else if (t >= mc2) {
      tref = mb2;
      asref = asmb;
      b0 = beta0(4) / (2. * M_PI);
    } else {
      tref = mc2;
      asref = asmc;
      b0 = beta0(3) / (2. * M_PI);
    }
    return 1. / (1. / asref + b0 * log(t / tref));
  }

  double as2(double t) const {
    /**
     * @brief calculate the strong coupling constant at 2 loops (next-to-leading
     * order)
     *
     * @param t the scale
     * @return the strong coupling constant
     */

    double tref, asref, b0, b1, w;

    if (t >= mb2) {
      tref = mz2;
      asref = asmz;
      b0 = beta0(5) / (2. * M_PI);
      b1 = beta1(5) / pow(2. * M_PI, 2);
    } else if (t >= mc2) {
      tref = mb2;
      asref = asmb;
      b0 = beta0(4) / (2. * M_PI);
      b1 = beta1(4) / pow(2. * M_PI, 2);
    } else {
      tref = mc2;
      asref = asmc;
      b0 = beta0(3) / (2. * M_PI);
      b1 = beta1(3) / pow(2. * M_PI, 2);
    }
    w = 1. + b0 * asref * log(t / tref);
    return asref / w * (1. - b1 / b0 * asref * log(w) / w);
  }

  double k_cmw(double t) const {
    /**
     * @brief calculate the CMW scheme conversion factor
     *
     * @param nf the number of flavours
     * @return the CMW scheme conversion factor
     */

    // Do threshold matching for the number of flavours
    int nf;
    if (t >= mb2) {
      nf = 5;
    } else if (t >= mc2) {
      nf = 4;
    } else {
      nf = 3;
    }

    // Do the numerical guards too
    nf = (fabs(t - mb2) < 1e-6 * mb2) ? 5 : nf;
    nf = (fabs(t - mc2) < 1e-6 * mc2) ? 4 : nf;

    return (67. / 18. - M_PI * M_PI / 6.) * k_ca - 10. / 9. * k_tr * nf;
  }

  double operator()(double t) const {
    /**
     * @brief wrapper/call operator for the strong coupling constant. This
     * function will calculate the strong coupling constant at the given scale
     * and at the number of loops specified in the constructor.
     *
     * @param t the scale
     * @return the strong coupling constant
     */

    // Freeze the scale to 1 GeV^2 if lower
    if (t < 1.) {
      t = 1.;
    }

    // Snap to threshold to avoid floating-point jitter at boundaries
    t = (fabs(t - mb2) < 1e-6 * mb2) ? mb2 : t;
    t = (fabs(t - mc2) < 1e-6 * mc2) ? mc2 : t;

    // First, calculate alpha_s
    double as_value;
    switch (n_loops) {
      case 0:
        as_value = asmz;  // No running coupling, for fixas tests
        break;
      case 1:
        as_value = as1(t);
        break;
      case 2:
        as_value = as2(t);
        break;
      default:
        as_value = as2(t);  // Default to 2-loop calculation
        break;
    }

    // Next, convert to CMW scheme if needed
    if (use_cmw) {
      double k_cmw_val = k_cmw(t);
      as_value *= 1. + k_cmw_val * as_value / (2. * M_PI);
    }
    return as_value;
  }
};

#endif  // qcd_h_