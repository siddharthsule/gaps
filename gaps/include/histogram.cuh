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

  __host__ __device__ Bin1D(double xmin = 0., double xmax = 0.)
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

class Bin2D {
 public:
  double xmin, xmax, ymin, ymax, w, w2, wx, wx2, wy, wy2, wxy, n;

 public:
  __host__ __device__ Bin2D(double xmin = 0., double xmax = 0.,
                            double ymin = 0., double ymax = 0.)
      : xmin(xmin),
        xmax(xmax),
        ymin(ymin),
        ymax(ymax),
        w(0.0),
        w2(0.0),
        wx(0.0),
        wx2(0.0),
        wy(0.0),
        wy2(0.0),
        wxy(0.0),
        n(0.0) {}

  __host__ __device__ double WidthX() const { return xmax - xmin; }
  __host__ __device__ double WidthY() const { return ymax - ymin; }

  __device__ void AtomicFill(double x, double y, double weight) {
    atomicAdd(&w, weight);
    atomicAdd(&w2, weight * weight);
    atomicAdd(&wx, weight * x);
    atomicAdd(&wx2, weight * weight * x);
    atomicAdd(&wy, weight * y);
    atomicAdd(&wy2, weight * weight * y);
    atomicAdd(&wxy, weight * x * y);
    atomicAdd(&n, 1.0);
  }

  __host__ __device__ void ScaleW(double scale) {
    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
    wy *= scale;
    wy2 *= scale * scale;
    wxy *= scale * scale;
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

class Histo2D {
 public:
  Bin2D bins[nBins][nBins];
  Bin2D uflow;
  Bin2D oflow;
  Bin2D total;
  double scale;

 public:
  __host__ __device__ Histo2D(double xmin = 0.0, double xmax = 1.0,
                              double ymin = 0.0, double ymax = 1.0)
      : uflow(xmin - 100., xmin, ymin - 100., ymin),
        oflow(xmax, xmax + 100., ymax, ymax + 100.),
        total(xmin - 100., xmax + 100., ymin - 100., ymax + 100.), scale(1.) {
    double widthX = (xmax - xmin) / nBins;
    double widthY = (ymax - ymin) / nBins;
    for (int i = 0; i < nBins; ++i) {
      for (int j = 0; j < nBins; ++j) {
        double xlow = xmin + i * widthX;
        double xhigh = xlow + widthX;
        double ylow = ymin + j * widthY;
        double yhigh = ylow + widthY;
        bins[i][j] = Bin2D(xlow, xhigh, ylow, yhigh);
      }
    }
  }

  __device__ void Fill(double x, double y, double w) {
    // Find the bin for the x-coordinate
    int lx = 0;
    int rx = nBins - 1;
    int cx = (lx + rx) / 2;
    double ax = bins[cx][0].xmin;

    while (rx - lx > 1) {
      if (x < ax) {
        rx = cx;
      } else {
        lx = cx;
      }
      cx = (lx + rx) / 2;
      ax = bins[cx][0].xmin;
    }

    // Find the bin for the y-coordinate
    int ly = 0;
    int ry = nBins - 1;
    int cy = (ly + ry) / 2;
    double ay = bins[0][cy].ymin;

    while (ry - ly > 1) {
      if (y < ay) {
        ry = cy;
      } else {
        ly = cy;
      }
      cy = (ly + ry) / 2;
      ay = bins[0][cy].ymin;
    }

    // Fill the appropriate bin
    if (x > bins[rx][0].xmin && y > bins[0][ry].ymin) {
      if (x > bins[rx][0].xmax || y > bins[0][ry].ymax) {
        oflow.AtomicFill(x, y, w);
      } else {
        bins[rx][ry].AtomicFill(x, y, w);
      }
    } else if (x < bins[lx][0].xmin || y < bins[0][ly].ymin) {
      uflow.AtomicFill(x, y, w);
    } else {
      bins[lx][ly].AtomicFill(x, y, w);
    }

    total.AtomicFill(x, y, w);
  }

  void ScaleW(double scale) {
    for (auto& binRow : bins) {
      for (auto& bin : binRow) {
        bin.ScaleW(scale);
      }
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

// Overload for Histo2D
std::string ToString(Histo2D h, std::string name);
void Write(Histo2D h, std::string name, const std::string& filename);

#endif  // HISTOGRAM_CUH_