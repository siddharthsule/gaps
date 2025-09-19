#pragma once
#ifndef LHAPDF_CuPDF_H
#define LHAPDF_CuPDF_H

#include <cuda.h>
#include <curand.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>

#include "CuInterpolator.h"
#include "LHAPDF/AlphaS.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/LHAPDF.h"

#define CUDA_LHAPDF

__device__ inline double ddlogq_forward(size_t i, double *alphas,
                                        double *logq2s) {
  return (alphas[i + 1] - alphas[i]) / (logq2s[i + 1] - logq2s[i]);
}

__device__ inline double ddlogq_backward(size_t i, double *alphas,
                                         double *logq2s) {
  return (alphas[i] - alphas[i - 1]) / (logq2s[i] - logq2s[i - 1]);
}

__device__ inline double ddlogq_central(size_t i, double *alphas,
                                        double *logq2s) {
  return 0.5 * (ddlogq_forward(i, alphas, logq2s) +
                ddlogq_backward(i, alphas, logq2s));
}

__global__ inline void _alphas(double *grid, double *values, double *q2s,
                               double *ret, size_t shape, int offset) {
  const int tid = (threadIdx.x + blockDim.x * blockIdx.x) + offset;
  double logq2 = log(q2s[tid]);
  if ((logq2 < grid[0]) || (logq2 > grid[shape - 1])) {
    ret[tid] = 0;
    return;
  }
  size_t iq2 = _binarySearch(logq2, grid, shape);
  double didlogq2, di1dlogq2;
  if (iq2 == 0 || (grid[iq2] == grid[iq2 - 1])) {
    didlogq2 = ddlogq_forward(iq2, values, grid);
    di1dlogq2 = ddlogq_central(iq2 + 1, values, grid);
  } else if (iq2 == shape - 2 || (grid[iq2 + 1] == grid[iq2 + 2])) {
    didlogq2 = ddlogq_central(iq2, values, grid);
    di1dlogq2 = ddlogq_backward(iq2 + 1, values, grid);
  } else {
    didlogq2 = ddlogq_central(iq2, values, grid);
    di1dlogq2 = ddlogq_central(iq2 + 1, values, grid);
  }

  // Calculate alpha_s
  const double dlogq2 = grid[iq2 + 1] - grid[iq2];
  const double tlogq2 = (logq2 - grid[iq2]) / dlogq2;
  const double v = _interpolateCubic(tlogq2, values[iq2], didlogq2 * dlogq2,
                                     values[iq2 + 1], di1dlogq2 * dlogq2);
  ret[tid] = v;
}

__global__ inline void _xfxQ2(double *grid, double *xknots, double *q2knots,
                              double *ret, size_t *shape, double *x, double *q2,
                              int pid, int *lookup, double *coeffs,
                              int offset) {
  const int tid = (threadIdx.x + blockDim.x * blockIdx.x) + offset;

  double logx = log(x[tid]);
  double logq2 = log(q2[tid]);
  int id = lookup[pid + 6];
  if (((logx < xknots[0]) or (logx > xknots[shape[0] - 1])) or
      ((logq2 < q2knots[0]) or (logq2 > q2knots[shape[1] - 1]))) {
    ret[tid] = 0;
    return;
  }
  size_t ix = _binarySearch(logx, xknots, shape[0]);
  size_t iq2 = _binarySearch(logq2, q2knots, shape[1]);

  // should probably sync somewhere here
  shared_data shared;
  shared.ix = ix;
  shared.iq2 = iq2;
  shared.grid = grid;
  shared.xknots = xknots;
  shared.q2knots = q2knots;
  shared.coeffs = coeffs;
  shared.shape = shape;
  shared.upper =
      ((iq2 + 1 == shape[1] - 1) || (q2knots[iq2 + 1] == q2knots[iq2 + 2]));
  shared.lower = ((iq2 == 0) || (q2knots[iq2] == q2knots[iq2 - 1]));

  ret[tid] = (id == -1) ? 0 : _interpolate(logx, logq2, id, shared);
}

class CuPDF {
 public:
  CuPDF(const std::string _setname, int _member)
      : setname(_setname), member(_member) {
    // init LHAPDF
    const std::string ipol = "logcubic";
    LHAPDF::GridPDF gridpdf(setname, member);
    gridpdf.setInterpolator(ipol);

    std::vector<size_t> &shape = gridpdf.Data().setShape();

    cudaMalloc((void **)&d_shape, sizeof(size_t) * 3);
    cudaMemcpy(d_shape, shape.data(), sizeof(size_t) * 3,
               cudaMemcpyHostToDevice);

    int n_grid_values = shape[0] * shape[1] * shape[2];
    int n_coeffs = (shape[0] - 1) * shape[1] * shape[2] * 4;
    cudaMalloc((void **)&d_xknots, sizeof(double) * shape[0]);
    cudaMalloc((void **)&d_q2knots, sizeof(double) * shape[1]);
    cudaMalloc((void **)&d_grid, sizeof(double) * n_grid_values);
    cudaMalloc((void **)&d_coeffs, sizeof(double) * n_coeffs);
    cudaMemcpy(d_xknots, gridpdf.Data().logxs().data(),
               sizeof(double) * shape[0], cudaMemcpyHostToDevice);
    cudaMemcpy(d_q2knots, gridpdf.Data().logq2s().data(),
               sizeof(double) * shape[1], cudaMemcpyHostToDevice);
    cudaMemcpy(d_grid, gridpdf.Data().setGrid().data(),
               sizeof(double) * n_grid_values, cudaMemcpyHostToDevice);
    cudaMemcpy(d_coeffs, gridpdf.Data().setCoeffs().data(),
               sizeof(double) * n_coeffs, cudaMemcpyHostToDevice);

    cudaMalloc((void **)&d_pids, sizeof(int) * shape[2]);
    cudaMalloc((void **)&d_lookup, sizeof(int) * 14);
    cudaMemcpy(d_pids, gridpdf.Data().setPids().data(), sizeof(int) * shape[2],
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_lookup, gridpdf.Data().setLookup().data(), sizeof(int) * 14,
               cudaMemcpyHostToDevice);

    // Alphas
    // TODO: Should have a sanity check here
    //  e.g. check if the alphas is actually a Ipol alphas
    const LHAPDF::AlphaS_Ipol &grid_alphas =
        dynamic_cast<LHAPDF::AlphaS_Ipol &>(gridpdf.alphaS());
    auto qs = grid_alphas.getQ2Values();
    // TODO: there must be a better solution for this..
    for (int i{0}; i < qs.size(); ++i) {
      qs[i] = log(qs[i]);
    }
    const auto as = grid_alphas.getAlphaSValues();

    // cuda memeory for alphas
    alphas_shape = qs.size();
    cudaMalloc((void **)&d_q2knots_as, sizeof(double) * qs.size());
    cudaMalloc((void **)&d_alphas, sizeof(double) * as.size());
    cudaMemcpy(d_q2knots_as, qs.data(), sizeof(double) * qs.size(),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_alphas, as.data(), sizeof(double) * as.size(),
               cudaMemcpyHostToDevice);
  }

  // Functions for computing the pdf for an array of memory
  void xfxQ2(double *xs, double *q2s, double *ret, int pid, int n_events) {
    // TODO implement all kinds of sanity checks
    int id = pid;
    if (pid == 21) id = 0;
    if (pid == 22) id = 14;
    if (n_events >= 256) {
      _xfxQ2<<<n_events / 256, 256>>>(d_grid, d_xknots, d_q2knots, ret, d_shape,
                                      xs, q2s, id, d_lookup, d_coeffs, 0);
    }
    if (n_events % 256 != 0) {
      const int offset{256 * (n_events / 256)};
      _xfxQ2<<<1, n_events % 256>>>(d_grid, d_xknots, d_q2knots, ret, d_shape,
                                    xs, q2s, id, d_lookup, d_coeffs, offset);
    }
  }

  void alphasQ2(double *q2s, double *ret, int n_events) {
    if (n_events >= 256) {
      _alphas<<<n_events / 256, 256>>>(d_q2knots_as, d_alphas, q2s, ret,
                                       alphas_shape, 0);
    }

    if (n_events % 256 != 0) {
      const int offset{256 * (n_events / 256)};
      _alphas<<<1, n_events % 256>>>(d_q2knots_as, d_alphas, q2s, ret,
                                     alphas_shape, offset);
    }
  }

 private:
  int blocksize, nblocks, member;
  std::string setname;

  // cuda memory for the PDF
  size_t *d_shape;
  double *d_xknots, *d_q2knots, *d_grid, *d_coeffs;
  int *d_pids, *d_lookup;

  // cuda memeory for alphas
  size_t alphas_shape;
  double *d_q2knots_as, *d_alphas;
};

#endif