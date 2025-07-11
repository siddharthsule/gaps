#ifndef histogram_h_
#define histogram_h_

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "base.h"

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

  bin1d(double xmin = 0., double xmax = 0.)
      : xmin(xmin), xmax(xmax), w(0.), w2(0.), wx(0.), wx2(0.), n(0.) {}

  // ---------------------------------------------------------------------------
  // member functions

  double width() const {
    /**
     * @brief get the width of the bin
     *
     * @return double: the width of the bin
     */

    return xmax - xmin;
  }

  void fill(double x, double weight) {
    /**
     * @brief fill the bin with a weight
     *
     * @param x the value to fill the bin with
     * @param weight the weight to fill the bin with
     */

    this->w += weight;
    w2 += weight * weight;
    wx += weight * x;
    wx2 += weight * weight * x;
    n += 1.;
  }

  void scale_w(double scale) {
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

  std::string format(const std::string& tag) const {
    /**
     * @brief format the bin data for the yoda file
     *
     * @param tag the tag for the bin
     */

    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << tag << "\t" << tag << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  std::string to_string() const {
    /**
     * @brief convert the bin data to a string
     */

    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << xmin << "\t" << xmax << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << static_cast<int>(n);
    return ss.str();
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

  std::string name;
  bin1d bins[n_bins];
  bin1d uflow;
  bin1d oflow;
  bin1d total;
  double scale;

  // ---------------------------------------------------------------------------
  // constructor

  histo1d(double xmin = 0., double xmax = 1., const std::string& name = "hst")
      : name(name),
        uflow(xmin - 100., xmin),
        oflow(xmax, xmax + 100.),
        total(xmin - 100., xmax + 100.),
        scale(1.) {
    double width = (xmax - xmin) / n_bins;
    for (int i = 0; i < n_bins; ++i) {
      double xlow = xmin + i * width;
      double xhigh = xlow + width;

      // Guards to prevent discontinuities in the histogram
      xlow = std::fabs(xlow) < 1e-12 ? 0. : xlow;
      xhigh = std::fabs(xhigh) < 1e-12 ? 0. : xhigh;

      // initialise the bins
      bins[i] = bin1d(xlow, xhigh);
    }
  }

  // ---------------------------------------------------------------------------
  // member functions

  void fill(double x, double w) {
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

    if (x > bins[r].xmin) {
      if (x > bins[r].xmax) {
        oflow.fill(x, w);
      } else {
        bins[r].fill(x, w);
      }
    } else if (x < bins[l].xmin) {
      uflow.fill(x, w);
    } else {
      bins[l].fill(x, w);
    }

    total.fill(x, w);
  }

  void scale_w(double scale) {
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

  std::string to_string() const {
    /**
     * @brief convert the histogram data to a string
     */

    std::stringstream ss;
    ss << "BEGIN YODA_HISTO1D " << name << "\n\n";
    ss << "Path=" << name << "\n\n";
    ss << "ScaledBy=" << scale << "\n";
    ss << "Title=\nType=Histo1D\n";
    ss << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    ss << total.format("Total") << "\n";
    ss << uflow.format("Underflow") << "\n";
    ss << oflow.format("Overflow") << "\n";
    ss << "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    for (const auto& bin : bins) {
      ss << bin.to_string() << "\n";
    }
    ss << "END YODA_HISTO1D\n\n";
    return ss.str();
  }

  void write(const std::string& filename) const {
    /**
     * @brief write the histogram data to a yoda file
     *
     * @param filename the name of the file to write the histogram data to
     */

    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);
    file << to_string();
    file.close();
  }
};

#endif  // histogram_h_