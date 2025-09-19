#ifndef histogram_cuh_
#define histogram_cuh_

#include "base.cuh"

class bin1d {
  /**
   * @class bin1d
   * @brief 1D bin for histogramming
   *
   * This class is used to store the data for a 1D bin in a histogram. It
   * contains the bin edges, the sum of weights, the sum of weights squared, the
   * sum of weights times the bin centre, the sum of weights squared times the
   * bin centre, and the number of entries.  the data is binned and then stored
   * as a yoda file[1]
   *
   * yoda: https://yoda.hepforge.org/
   */

 public:
  // ---------------------------------------------------------------------------
  // member variables

  double xmin, xmax, w, w2, wx, wx2, n;

  // ---------------------------------------------------------------------------
  // constructor

  __host__ __device__ bin1d(double xmin = 0., double xmax = 0.)
      : xmin(xmin), xmax(xmax), w(0.), w2(0.), wx(0.), wx2(0.), n(0.) {}

  // ---------------------------------------------------------------------------
  // member functions

  __host__ __device__ double width() const {
    /**
     * @brief get the width of the bin
     *
     * @return double: the width of the bin
     */

    return xmax - xmin;
  }

  __device__ void atomic_fill(double x, double weight) {
    /**
     * @brief fill the bin with a weight
     *
     * @param x the value to fill the bin with
     * @param weight the weight to fill the bin with
     */

    atomicAdd(&w, weight);
    atomicAdd(&w2, weight * weight);
    atomicAdd(&wx, weight * x);
    atomicAdd(&wx2, weight * weight * x);
    atomicAdd(&n, 1.);
  }

  __host__ __device__ void scale_w(double scale) {
    /**
     * @brief scale the weight of the bin
     *
     * @param scale the scale factor
     */

    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
  }
};

// ---------------------------------------------------------------------------

class histo1d {
  /**
   * @class histo1d
   * @brief 1D histogram
   *
   * This class is used to store the data for a 1D histogram. It contains the
   * name of the histogram, the bins, the underflow bin, the overflow bin, the
   * total bin, and the scale factor. The data is binned and then stored as a
   * yoda file[1]
   *
   * yoda: https://yoda.hepforge.org/
   */

 public:
  // ---------------------------------------------------------------------------
  // member variables

  /**
   * name component of histo1d
   * -------------------------
   *
   * unfortunately, std::string is not a feature in cuda, so we have to provide
   * the name additinally whern writing to file. this is the only difference
   * between the GPU and CPU versions of the code.
   */

  char name[20];
  bin1d bins[max_bins];
  int n_bins;  // Actual number of bins used
  bin1d uflow;
  bin1d oflow;
  bin1d total;
  double scale;

  // ---------------------------------------------------------------------------
  // constructor

  __host__ __device__ histo1d(int n_bins = max_bins, double xmin = 0.,
                              double xmax = 1.,
                              const char* hist_name = "histo1d",
                              bool logspace = false)
      : n_bins(n_bins),
        uflow(xmin - 100., xmin),
        oflow(xmax, xmax + 100.),
        total(xmin - 100., xmax + 100.),
        scale(1.) {
    // copy the name to the name variable
    for (int i = 0; i < 19; ++i) {
      name[i] = hist_name[i];
      if (hist_name[i] == '\0') break;
    }
    name[19] = '\0';  // Ensure null-termination

    // Option 1: Linear spacing
    if (!logspace) {
      double width = (xmax - xmin) / n_bins;
      for (int i = 0; i < n_bins; ++i) {
        double xlow = xmin + i * width;
        double xhigh = xlow + width;

        // Guards to prevent discontinuities in the histogram
        xlow = abs(xlow) < 1e-12 ? 0. : xlow;
        xhigh = abs(xhigh) < 1e-12 ? 0. : xhigh;

        // initialize bin1d object on the device
        bins[i] = bin1d(xlow, xhigh);
      }
    }

    // Option 2: Logarithmic spacing
    else {
      double log_xmin = log10(xmin);
      double log_xmax = log10(xmax);
      double log_width = (log_xmax - log_xmin) / n_bins;

      for (int i = 0; i < n_bins; ++i) {
        double xlow = pow(10, log_xmin + static_cast<double>(i) * log_width);
        double xhigh =
            pow(10, log_xmin + static_cast<double>(i + 1) * log_width);

        // Guards to prevent discontinuities in the histogram
        xlow = fabs(xlow) < 1e-12 ? 0. : xlow;
        xhigh = fabs(xhigh) < 1e-12 ? 0. : xhigh;

        // initialize bin1d object on the device
        bins[i] = bin1d(xlow, xhigh);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // member functions

  __device__ void fill(double x, double w) {
    /**
     * @brief fill the histogram with a weight
     *
     * @param x the value to fill the histogram with
     * @param w the weight to fill the histogram with
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

    // atomic binning! each event is binned simultaneously here
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
    /**
     * @brief scale the weight of the histogram
     *
     * @param scale the scale factor
     */

    for (int i = 0; i < n_bins; ++i) {
      bins[i].scale_w(scale);
    }
    uflow.scale_w(scale);
    oflow.scale_w(scale);
    total.scale_w(scale);
    this->scale *= scale;
  }
};

// writing is done outside of the class in cuda implementation
std::string to_string(histo1d h);
void write(histo1d h, const std::string& filename);

#endif  // histogram_cuh_