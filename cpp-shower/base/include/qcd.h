#ifndef QCD_H_
#define QCD_H_

// Base Class, with all the important definitions
#include "base.h"

const double kNC = 3.0;
const double kTR = 0.5;
const double kCA = kNC;
const double kCF = (kNC * kNC - 1.0) / (2.0 * kNC);

class AlphaS {
 private:
  int order;
  double mc2, mb2, mz2, asmz, asmb, asmc;

 public:
  // Constructor
  AlphaS(double mz, double asmz, int order = 1, double mb = 4.75,
         double mc = 1.27)
      : order(order),
        mc2(mc * mc),
        mb2(mb * mb),
        mz2(mz * mz),
        asmz(asmz),
        asmb((*this)(mb2)),
        asmc((*this)(mc2)) {}

  // Beta functions
  double Beta0(int nf) const { return (11. / 6. * kCA) - (2. / 3. * kTR * nf); }

  double Beta1(int nf) const {
    return (17. / 6. * kCA * kCA) - ((5. / 3. * kCA + kCF) * kTR * nf);
  }

  // Alpha_s at order 0 and 1 (One-Loop and Two-Loop)
  double As0(double t) const {
    double tref, asref, b0;
    if (t >= mb2) {
      tref = mz2;
      asref = asmz;
      b0 = Beta0(5) / (2. * M_PI);
    } else if (t >= mc2) {
      tref = mb2;
      asref = asmb;
      b0 = Beta0(4) / (2. * M_PI);
    } else {
      tref = mc2;
      asref = asmc;
      b0 = Beta0(3) / (2. * M_PI);
    }
    return 1. / (1. / asref + b0 * log(t / tref));
  }

  double As1(double t) const {
    double tref, asref, b0, b1, w;
    if (t >= mb2) {
      tref = mz2;
      asref = asmz;
      b0 = Beta0(5) / (2. * M_PI);
      b1 = Beta1(5) / pow(2. * M_PI, 2);
    } else if (t >= mc2) {
      tref = mb2;
      asref = asmb;
      b0 = Beta0(4) / (2. * M_PI);
      b1 = Beta1(4) / pow(2. * M_PI, 2);
    } else {
      tref = mc2;
      asref = asmc;
      b0 = Beta0(3) / (2. * M_PI);
      b1 = Beta1(3) / pow(2. * M_PI, 2);
    }
    w = 1. + b0 * asref * log(t / tref);
    return asref / w * (1. - b1 / b0 * asref * log(w) / w);
  }

  // Call operator to calculate alpha_s
  double operator()(double t) {
    if (order == 0) {
      return As0(t);
    } else {
      return As1(t);
    }
  }
};

#endif  // QCD_H_