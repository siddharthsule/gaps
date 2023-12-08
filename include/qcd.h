#ifndef QCD_H_
#define QCD_H_

#include <cmath>

class AlphaS {
 public:
  // Constructor
  AlphaS(double mz, double asmz, int order = 1, double mb = 4.75, double mc = 1.27);

  // Beta functions
  double Beta0(int nf);
  double Beta1(int nf);

  // Alpha_s at order 0 and 1
  double As0(double t);
  double As1(double t);

  // Call operator to calculate alpha_s
  double operator()(double t);

 private:
  int order;
  double mc2, mb2, mz2, asmz, asmb, asmc;
};

const double kNC = 3.0;
const double kTR = 0.5;
const double kCA = kNC;
const double kCF = (kNC * kNC - 1.0) / (2.0 * kNC);

#endif  // QCD_H_