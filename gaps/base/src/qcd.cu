#include "qcd.cuh"

// constructor
__device__ alpha_s::alpha_s(double mz, double asmz, int order, double mb,
                            double mc)
    : order(order),
      mc2(mc * mc),
      mb2(mb * mb),
      mz2(mz * mz),
      asmz(asmz),
      asmb((*this)(mb2)),
      asmc((*this)(mc2)) {}

// setup
__device__ void alpha_s::setup(double mz, double asmz, int order, double mb,
                               double mc) {
  this->order = order;
  this->mc2 = mc * mc;
  this->mb2 = mb * mb;
  this->mz2 = mz * mz;
  this->asmz = asmz;
  this->asmb = (*this)(mb2);
  this->asmc = (*this)(mc2);
}

// beta and alpha s functions
__device__ double alpha_s::beta0(int nf) {
  return (11. / 6. * k_ca) - (2. / 3. * k_tr * nf);
}

__device__ double alpha_s::beta1(int nf) {
  return (17. / 6. * k_ca * k_ca) - ((5. / 3. * k_ca + k_cf) * k_tr * nf);
}

// alpha_s at order 0 and 1 (one-loop and two-loop)
__device__ double alpha_s::as0(double t) {
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

__device__ double alpha_s::as1(double t) {
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

__device__ double alpha_s::operator()(double t) {
  if (order == 0) {
    return as0(t);
  } else {
    return as1(t);
  }
}

// set up kernel on the device
__global__ void as_setup_kernel(alpha_s *as, double mz, double asmz, int order,
                                double mb, double mc) {
  as->setup(mz, asmz, order, mb, mc);
}

// calculate alpha_s on the device for one input
__global__ void as_value(alpha_s *as, double *asval, double t) {
  int idx = block_idx.x * block_dim.x + thread_idx.x;

  if (idx >= 1) return;

  asval[idx] = (*as)(t);
  printf("as_val: %f\n", (*as)(t));
}

// calculate alpha_s on the device for many inputs
// exclusively used for parton shower veto algorithm
__global__ void as_shower_kernel(alpha_s *as, event *events, double *asval,
                                 int n) {
  int idx = block_idx.x * block_dim.x + thread_idx.x;

  if (idx >= n) return;
  event &ev = events[idx];

  asval[idx] = (*as)(ev.get_shower_t());
}
