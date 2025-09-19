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
  double mc2, mb2, mz2, asmz, asmb, asmc;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  alpha_s(double mz, double asmz, int n_loops = 2, double mb = 4.75,
          double mc = 1.3)
      : n_loops(n_loops),
        mc2(mc * mc),
        mb2(mb * mb),
        mz2(mz * mz),
        asmz(asmz),
        asmb((*this)(mb2)),
        asmc((*this)(mc2)) {}

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

    // Threshold Matching is OFF - massless quarks!
    // if (t >= mb2) {
    tref = mz2;
    asref = asmz;
    b0 = beta0(5) / (2. * M_PI);
    // } else if (t >= mc2) {
    //   tref = mb2;
    //   asref = asmb;
    //   b0 = beta0(4) / (2. * M_PI);
    // } else {
    //   tref = mc2;
    //   asref = asmc;
    //   b0 = beta0(3) / (2. * M_PI);
    // }
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

    // Threshold Matching is OFF - massless quarks!
    // if (t >= mb2) {
    tref = mz2;
    asref = asmz;
    b0 = beta0(5) / (2. * M_PI);
    b1 = beta1(5) / pow(2. * M_PI, 2);
    // } else if (t >= mc2) {
    //   tref = mb2;
    //   asref = asmb;
    //   b0 = beta0(4) / (2. * M_PI);
    //   b1 = beta1(4) / pow(2. * M_PI, 2);
    // } else {
    //   tref = mc2;
    //   asref = asmc;
    //   b0 = beta0(3) / (2. * M_PI);
    //   b1 = beta1(3) / pow(2. * M_PI, 2);
    // }
    w = 1. + b0 * asref * log(t / tref);
    return asref / w * (1. - b1 / b0 * asref * log(w) / w);
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

    switch (n_loops) {
      case 0:
        return asmz;  // No running coupling, for fixas tests
      case 1:
        return as1(t);
      case 2:
        return as2(t);
      default:
        return as2(t);  // Default to 2-loop calculation
    }
  }
};

#endif  // qcd_h_