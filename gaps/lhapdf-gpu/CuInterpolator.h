#pragma once
#ifndef LHAPDF_CuInterpolator_H
#define LHAPDF_CuInterpolator_H

#include <cuda.h>
#include <curand.h>

#include "Cu1Dinterpolators.h"

struct shared_data {
  int ix, iq2;
  double *grid, *xknots, *q2knots, *coeffs;
  size_t* shape;
  bool upper, lower;
};
__device__ inline double xf(double* grid, size_t* shape, size_t ix, size_t iq2,
                            size_t pid) {
  return grid[ix * shape[2] * shape[1] + iq2 * shape[2] + pid];
}

__device__ inline int coeffIndex(int ix, int iq2, int id, size_t* shape) {
  return ix * shape[1] * shape[2] * 4 + iq2 * shape[2] * 4 + id * 4;
}

__device__ inline size_t _binarySearch(double x, double* knots, size_t _r) {
  int l = 0;
  int r = _r;
  while (l <= r) {
    size_t m = l + (r - l) / 2;
    if (knots[m] <= x && x <= knots[m + 1]) return m;
    if (knots[m] < x)
      l = m + 1;
    else
      r = m - 1;
  }
  return 0;
}

#ifdef LHAPDF_DEBUG
// logbilinear interpolation, mostly here for debugging
__device__ inline double _interpolate(double logx, double logq2, int id,
                                      shared_data shared) {
  int ix = shared.ix;
  int iq2 = shared.iq2;
  const double f_ql =
      _interpolateLinear(logx, shared.xknots[ix], shared.xknots[ix + 1],
                         xf(shared.grid, shared.shape, ix, iq2, id),
                         xf(shared.grid, shared.shape, ix + 1, iq2, id));
  const double f_qh =
      _interpolateLinear(logx, shared.xknots[ix], shared.xknots[ix + 1],
                         xf(shared.grid, shared.shape, ix, iq2 + 1, id),
                         xf(shared.grid, shared.shape, ix + 1, iq2 + 1, id));
  return _interpolateLinear(logq2, shared.q2knots[iq2], shared.q2knots[iq2 + 1],
                            f_ql, f_qh);
}
#else

__device__ inline double _interpolate(double logx, double logq2, int id,
                                      shared_data shared) {
  int ix = shared.ix;
  int iq2 = shared.iq2;

  double dlogx_1 = shared.xknots[ix + 1] - shared.xknots[ix];
  double tlogx = (logx - shared.xknots[ix]) / dlogx_1;
  double dlogq_0 = 1. / (shared.q2knots[iq2] - shared.q2knots[iq2 - 1]);
  double dlogq_1 = shared.q2knots[iq2 + 1] - shared.q2knots[iq2];
  double dlogq_2 = 1. / (shared.q2knots[iq2 + 2] - shared.q2knots[iq2 + 1]);
  double tlogq = (logq2 - shared.q2knots[iq2]) / dlogq_1;

  double vl = _interpolateCubic(tlogx, shared.coeffs,
                                coeffIndex(ix, iq2, id, shared.shape));
  double vh = _interpolateCubic(tlogx, shared.coeffs,
                                coeffIndex(ix, iq2 + 1, id, shared.shape));
  double vdl, vdh;
  if (shared.lower) {
    vdl = (vh - vl);
    double vhh = _interpolateCubic(tlogx, shared.coeffs,
                                   coeffIndex(ix, iq2 + 2, id, shared.shape));
    vdh = (vdl + (vhh - vh) * dlogq_1 * dlogq_2) * 0.5;
  } else if (shared.upper) {
    vdh = (vh - vl);
    double vll = _interpolateCubic(tlogx, shared.coeffs,
                                   coeffIndex(ix, iq2 - 1, id, shared.shape));
    vdl = (vdh + (vl - vll) * dlogq_1 * dlogq_0) * 0.5;
  } else {
    double vll = _interpolateCubic(tlogx, shared.coeffs,
                                   coeffIndex(ix, iq2 - 1, id, shared.shape));
    vdl = ((vh - vl) + (vl - vll) * dlogq_1 * dlogq_0) * 0.5;
    double vhh = _interpolateCubic(tlogx, shared.coeffs,
                                   coeffIndex(ix, iq2 + 2, id, shared.shape));
    vdh = ((vh - vl) + (vhh - vh) * dlogq_1 * dlogq_2) * 0.5;
  }
  return _interpolateCubic(tlogq, vl, vdl, vh, vdh);
}
#endif

#endif