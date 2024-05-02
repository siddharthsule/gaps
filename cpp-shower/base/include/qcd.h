#ifndef qcd_h_
#define qcd_h_

// base class, with all the important definitions
#include "base.h"

const double k_nc = 3.;
const double k_tr = 0.5;
const double k_ca = k_nc;
const double k_cf = (k_nc * k_nc - 1.) / (2. * k_nc);
class alpha_s {
 private:
  int order;
  double mc2, mb2, mz2, asmz, asmb, asmc;

 public:
  // constructor
  alpha_s(double mz, double asmz, int order = 1, double mb = 4.75,
          double mc = 1.27)
      : order(order),
        mc2(mc * mc),
        mb2(mb * mb),
        mz2(mz * mz),
        asmz(asmz),
        asmb((*this)(mb2)),
        asmc((*this)(mc2)) {}

  // beta functions
  double beta0(int nf) const {
    return (11. / 6. * k_ca) - (2. / 3. * k_tr * nf);
  }

  double beta1(int nf) const {
    return (17. / 6. * k_ca * k_ca) - ((5. / 3. * k_ca + k_cf) * k_tr * nf);
  }

  // alpha_s at order 0 and 1 (one-loop and two-loop)
  double as0(double t) const {
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

  double as1(double t) const {
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

  // call operator to calculate alpha_s
  double operator()(double t) {
    if (order == 0) {
      return as0(t);
    } else {
      return as1(t);
    }
  }
};

#endif  // qcd_h_