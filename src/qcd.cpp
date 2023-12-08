#include "qcd.h"

// Constructor
AlphaS::AlphaS(double mz, double asmz, int order, double mb, double mc) :
    order(order), mc2(mc * mc), mb2(mb * mb), mz2(mz * mz), asmz(asmz), asmb((*this)(mb2)), asmc((*this)(mc2)) {}

// Beta functions
double AlphaS::Beta0(int nf) {
  return (11. / 6. * kCA) - (2. / 3. * kTR * nf);
}

double AlphaS::Beta1(int nf) {
  return (17. / 6. * kCA * kCA) - ((5. / 3. * kCA + kCF) * kTR * nf);
}

// Alpha_s at order 0 and 1
double AlphaS::As0(double t) {
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

double AlphaS::As1(double t) {
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
double AlphaS::operator()(double t) {
  if (order == 0) {
    return As0(t);
  } else {
    return As1(t);
  }
}