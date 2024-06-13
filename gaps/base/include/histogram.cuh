#ifndef histogram_cuh_
#define histogram_cuh_

#include "base.cuh"

/**
 * binning and histogramming
 * -------------------------
 *
 * this file contains tools needed for binning and histogramming data. the data
 * is binned and then stored as a yoda file[1]
 *
 * yoda: https://yoda.hepforge.org/
 */

class bin1d {
 public:
  double xmin, xmax, w, w2, wx, wx2, n;

  __host__ __device__ bin1d(double xmin = 0., double xmax = 0.)
      : xmin(xmin), xmax(xmax), w(0.), w2(0.), wx(0.), wx2(0.), n(0.) {}

  __host__ __device__ double width() const { return xmax - xmin; }

  __device__ void atomic_fill(double x, double weight) {
    atomicAdd(&w, weight);
    atomicAdd(&w2, weight * weight);
    atomicAdd(&wx, weight * x);
    atomicAdd(&wx2, weight * weight * x);
    atomicAdd(&n, 1.);
  }

  __host__ __device__ void scale_w(double scale) {
    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
  }
};

class bin2d {
 public:
  double xmin, xmax, ymin, ymax, w, w2, wx, wx2, wy, wy2, wxy, n;

 public:
  __host__ __device__ bin2d(double xmin = 0., double xmax = 0.,
                            double ymin = 0., double ymax = 0.)
      : xmin(xmin),
        xmax(xmax),
        ymin(ymin),
        ymax(ymax),
        w(0.),
        w2(0.),
        wx(0.),
        wx2(0.),
        wy(0.),
        wy2(0.),
        wxy(0.),
        n(0.) {}

  __host__ __device__ double width_x() const { return xmax - xmin; }
  __host__ __device__ double width_y() const { return ymax - ymin; }

  __device__ void atomic_fill(double x, double y, double weight) {
    atomicAdd(&w, weight);
    atomicAdd(&w2, weight * weight);
    atomicAdd(&wx, weight * x);
    atomicAdd(&wx2, weight * weight * x);
    atomicAdd(&wy, weight * y);
    atomicAdd(&wy2, weight * weight * y);
    atomicAdd(&wxy, weight * x * y);
    atomicAdd(&n, 1.);
  }

  __host__ __device__ void scale_w(double scale) {
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
 * name component of histo1d
 * -------------------------
 *
 * unfortunately, std::string is not a feature in cuda, so we have to provide
 * the name additinally whern writing to file. this is the only difference
 * between the cuda and c++ versions of the code.
 */
// histo1d class
class histo1d {
 public:
  bin1d bins[n_bins];  // array of bin1d objects on the device
  bin1d uflow;
  bin1d oflow;
  bin1d total;
  double scale;
  // static constexpr int nbin = n_bins;  // set in base.cuh

 public:
  __host__ __device__ histo1d(double xmin = 0., double xmax = 1.)
      : uflow(xmin - 100., xmin),
        oflow(xmax, xmax + 100.),
        total(xmin - 100., xmax + 100.),
        scale(1.) {
    double width = (xmax - xmin) / n_bins;
    for (int i = 0; i < n_bins; ++i) {
      double xlow = xmin + i * width;
      double xhigh = xlow + width;
      bins[i] = bin1d(xlow, xhigh);  // initialize bin1d object on the device
    }
  }

  __device__ void fill(double x, double w) {
    /**
     * atomic binning! each event is binned simultaneously here
     */
    int l = 0;
    int r = n_bins - 1;
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
        oflow.atomic_fill(x, w);
      } else {
        bins[r].atomic_fill(x, w);
      }
    } else if (x < bins[l].xmin) {
      uflow.atomic_fill(x, w);
    } else {
      bins[l].atomic_fill(x, w);
    }

    total.atomic_fill(x, w);
  }

  __host__ __device__ void scale_w(double scale) {
    for (int i = 0; i < n_bins; ++i) {
      bins[i].scale_w(scale);
    }
    uflow.scale_w(scale);
    oflow.scale_w(scale);
    total.scale_w(scale);
    this->scale *= scale;
  }
};

class histo2d {
 public:
  bin2d bins[n_bins2d][n_bins2d];
  bin2d uflow;
  bin2d oflow;
  bin2d total;
  double scale;

 public:
  __host__ __device__ histo2d(double xmin = 0., double xmax = 1.,
                              double ymin = 0., double ymax = 1.)
      : uflow(xmin - 100., xmin, ymin - 100., ymin),
        oflow(xmax, xmax + 100., ymax, ymax + 100.),
        total(xmin - 100., xmax + 100., ymin - 100., ymax + 100.),
        scale(1.) {
    double width_x = (xmax - xmin) / n_bins2d;
    double width_y = (ymax - ymin) / n_bins2d;
    for (int i = 0; i < n_bins2d; ++i) {
      for (int j = 0; j < n_bins2d; ++j) {
        double xlow = xmin + i * width_x;
        double xhigh = xlow + width_x;
        double ylow = ymin + j * width_y;
        double yhigh = ylow + width_y;
        bins[i][j] = bin2d(xlow, xhigh, ylow, yhigh);
      }
    }
  }

  __device__ void fill(double x, double y, double w) {
    // find the bin for the x-coordinate
    int lx = 0;
    int rx = n_bins2d - 1;
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

    // find the bin for the y-coordinate
    int ly = 0;
    int ry = n_bins2d - 1;
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

    // fill the appropriate bin
    if (x > bins[rx][0].xmin && y > bins[0][ry].ymin) {
      if (x > bins[rx][0].xmax || y > bins[0][ry].ymax) {
        oflow.atomic_fill(x, y, w);
      } else {
        bins[rx][ry].atomic_fill(x, y, w);
      }
    } else if (x < bins[lx][0].xmin || y < bins[0][ly].ymin) {
      uflow.atomic_fill(x, y, w);
    } else {
      bins[lx][ly].atomic_fill(x, y, w);
    }

    total.atomic_fill(x, y, w);
  }

  void scale_w(double scale) {
    for (auto& bin_row : bins) {
      for (auto& bin : bin_row) {
        bin.scale_w(scale);
      }
    }
    uflow.scale_w(scale);
    oflow.scale_w(scale);
    total.scale_w(scale);
    this->scale *= scale;
  }
};

// writing is done outside of the class in cuda implementation
std::string to_string(histo1d h, std::string name);
void write(histo1d h, std::string name, const std::string& filename);

// overload for histo2d
std::string to_string(histo2d h, std::string name);
void write(histo2d h, std::string name, const std::string& filename);

#endif  // histogram_cuh_