#pragma once
#ifndef LHAPDF_Cu1Dinterpolators_H
#define LHAPDF_Cu1Dinterpolators_H

__device__ inline double _interpolateLinear(const double x, const double xl,
                                            const double xh, const double yl,
                                            const double yh) {
  return yl + (x - xl) / (xh - xl) * (yh - yl);
}

__device__ inline double _interpolateCubic(double x, double* coeffs,
                                           size_t index) {
  const double x2 = x * x;
  const double x3 = x2 * x;
  return coeffs[index + 0] * x3 + coeffs[index + 1] * x2 +
         coeffs[index + 2] * x + coeffs[index + 3];
}

__device__ inline double _interpolateCubic(double T, double VL, double VDL,
                                           double VH, double VDH) {
  // Pre-calculate powers of T
  const double t2 = T * T;
  const double t3 = t2 * T;

  // Calculate left point
  const double p0 = (2 * t3 - 3 * t2 + 1) * VL;
  const double m0 = (t3 - 2 * t2 + T) * VDL;

  // Calculate right point
  const double p1 = (-2 * t3 + 3 * t2) * VH;
  const double m1 = (t3 - t2) * VDH;

  return p0 + m0 + p1 + m1;
}

#endif