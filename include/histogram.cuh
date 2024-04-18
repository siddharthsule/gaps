#ifndef HISTOGRAM_CUH_
#define HISTOGRAM_CUH_

#include "base.cuh"

/**
 * Binning and Histogramming
 * -------------------------
 *
 * This file contains tools needed for binning and histogramming data. The data
 * is binned and then stored as a Yoda file[1]
 *
 * Yoda: https://yoda.hepforge.org/
 */

class Bin1D {
 public:
  double xmin, xmax, w, w2, wx, wx2, n;

  __host__ __device__ Bin1D(double xmin = 0, double xmax = 0)
      : xmin(xmin), xmax(xmax), w(0.0), w2(0.0), wx(0.0), wx2(0.0), n(0.0) {}

  __host__ __device__ double Width() const { return xmax - xmin; }

  __device__ void AtomicFill(double x, double weight) {
    atomicAdd(&w, weight);
    atomicAdd(&w2, weight * weight);
    atomicAdd(&wx, weight * x);
    atomicAdd(&wx2, weight * weight * x);
    atomicAdd(&n, 1.0);
  }

  __host__ __device__ void ScaleW(double scale) {
    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
  }
};

/**
 * Name component of Histo1D
 * -------------------------
 *
 * Unfortunately, std::string is not a feature in CUDA, so we have to provide
 * the name additinally whern writing to file. This is the only difference
 * between the CUDA and C++ versions of the code.
 */
// Histo1D class
class Histo1D {
 public:
  Bin1D bins[nBins];  // Array of Bin1D objects on the device
  Bin1D uflow;
  Bin1D oflow;
  Bin1D total;
  double scale;
  // static constexpr int nbin = nBins;  // SET IN BASE.CUH

 public:
  __host__ __device__ Histo1D(double xmin = 0.0, double xmax = 1.0)
      : uflow(xmin - 100., xmin),
        oflow(xmax, xmax + 100.),
        total(xmin - 100., xmax + 100.),
        scale(1.) {
    double width = (xmax - xmin) / nBins;
    for (int i = 0; i < nBins; ++i) {
      double xlow = xmin + i * width;
      double xhigh = xlow + width;
      bins[i] = Bin1D(xlow, xhigh);  // Initialize Bin1D object on the device
    }
  }

  // Atomic Binning! Each event is binned simultaneously here
  __device__ void Fill(double x, double w) {
    int l = 0;
    int r = nBins - 1;
    int c = (l + r) / 2;
    double a = bins[c].xmin;

    while (r - l > 1) {
      if (x < a) {
        r = c;
      } else {
        l = c;
      }
      c = (l + r) / 2;
      a = bins[c].xmin;
    }

    if (x > bins[r].xmin) {
      if (x > bins[r].xmax) {
        oflow.AtomicFill(x, w);
      } else {
        bins[r].AtomicFill(x, w);
      }
    } else if (x < bins[l].xmin) {
      uflow.AtomicFill(x, w);
    } else {
      bins[l].AtomicFill(x, w);
    }

    total.AtomicFill(x, w);
  }

  __host__ __device__ void ScaleW(double scale) {
    for (int i = 0; i < nBins; ++i) {
      bins[i].ScaleW(scale);
    }
    uflow.ScaleW(scale);
    oflow.ScaleW(scale);
    total.ScaleW(scale);
    this->scale *= scale;
  }
};

// Writing is done outside of the class in CUDA implementation
std::string ToString(Histo1D h, std::string name);
void Write(Histo1D h, std::string name, const std::string& filename);

#endif  // HISTOGRAM_CUH_