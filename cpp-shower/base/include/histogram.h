#ifndef histogram_h_
#define histogram_h_

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "base.h"

/**
 * binning and histogramming
 * -------------------------
 *
 * this file contains tools needed for binning and histogramming data. the data
 * is binned and then stored as a yoda file[1]
 *
 * yoda: https://yoda.hepforge.org/
 */

// bin1d class
class bin1d {
 public:
  double xmin, xmax, w, w2, wx, wx2, n;

 public:
  bin1d(double xmin, double xmax)
      : xmin(xmin), xmax(xmax), w(0.), w2(0.), wx(0.), wx2(0.), n(0.) {}

  std::string format(const std::string& tag) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << tag << "\t" << tag << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  std::string to_string() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << xmin << "\t" << xmax << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  double width() const { return xmax - xmin; }

  void fill(double x, double weight) {
    this->w += weight;
    w2 += weight * weight;
    wx += weight * x;
    wx2 += weight * weight * x;
    n += 1.;
  }

  void scale_w(double scale) {
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
  bin2d(double xmin, double xmax, double ymin, double ymax)
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

  std::string format(const std::string& tag) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << tag << "\t" << tag << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << wy << "\t" << wy2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  std::string to_string() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << xmin << "\t" << xmax << "\t" << ymin << "\t" << ymax << "\t" << w
       << "\t" << w2 << "\t" << wx << "\t" << wx2 << "\t" << wy << "\t" << wy2
       << "\t" << wxy << "\t" << static_cast<int>(n);
    return ss.str();
  }

  double width_x() const { return xmax - xmin; }
  double width_y() const { return ymax - ymin; }

  void fill(double x, double y, double weight) {
    this->w += weight;
    w2 += weight * weight;
    wx += weight * x;
    wx2 += weight * weight * x;
    wy += weight * y;
    wy2 += weight * weight * y;
    wxy += weight * x * y;
    n += 1.;
  }

  void scale_w(double scale) {
    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
    wy *= scale;
    wy2 *= scale * scale;
    wxy *= scale * scale;
  }
};

// histo1d class
class histo1d {
 public:
  std::string name;
  std::vector<bin1d> bins;
  bin1d uflow;
  bin1d oflow;
  bin1d total;
  double scale;

 public:
  // constructor for histo1d
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
      bins.push_back(bin1d(xlow, xhigh));
    }
  }

  std::string to_string() const {
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

  void fill(double x, double w) {
    int l = 0;
    int r = bins.size() - 1;
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
    for (auto& bin : bins) {
      bin.scale_w(scale);
    }
    uflow.scale_w(scale);
    oflow.scale_w(scale);
    total.scale_w(scale);
    this->scale *= scale;
  }

  void write(const std::string& filename) const {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);
    file << to_string();
    file.close();
  }
};

// histo2d class
class histo2d {
 public:
  std::string name;
  std::vector<std::vector<bin2d>> bins;
  bin2d uflow;
  bin2d oflow;
  bin2d total;
  double scale;

 public:
  // constructor for histo2d
  histo2d(double xmin = 0., double xmax = 1., double ymin = 0.,
          double ymax = 1., const std::string& name = "hst")
      : name(name),
        uflow(xmin - 100., xmin, ymin - 100., ymin),
        oflow(xmax, xmax + 100., ymax, ymax + 100.),
        total(xmin - 100., xmax + 100., ymin - 100., ymax + 100.),
        scale(1.) {
    double xwidth = (xmax - xmin) / n_bins2d;
    double ywidth = (ymax - ymin) / n_bins2d;
    for (int i = 0; i < n_bins2d; ++i) {
      double xlow = xmin + i * xwidth;
      double xhigh = xlow + xwidth;
      std::vector<bin2d> bin_row;
      for (int j = 0; j < n_bins2d; ++j) {
        double ylow = ymin + j * ywidth;
        double yhigh = ylow + ywidth;
        bin_row.push_back(bin2d(xlow, xhigh, ylow, yhigh));
      }
      bins.push_back(bin_row);
    }
  }

  std::string to_string() const {
    std::stringstream ss;
    ss << "BEGIN YODA_HISTO2D " << name << "\n\n";
    ss << "Path=" << name << "\n\n";
    ss << "ScaledBy=" << scale << "\n";
    ss << "Title=\nType=Histo2D\n";
    ss << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2\tnumEntries\n";
    ss << total.format("Total") << "\n";
    ss << "# "
          "xlow\txhigh\tylow\tyhigh\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2"
          "\tnumEntries\n";
    for (const auto& binRow : bins) {
      for (const auto& bin : binRow) {
        ss << bin.to_string() << "\n";
      }
    }
    ss << "END YODA_HISTO2D\n\n";
    return ss.str();
  }

  void fill(double x, double y, double w) {
    // find the bin for the x-coordinate
    int lx = 0;
    int rx = bins.size() - 1;
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
    int ry = bins[0].size() - 1;
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
        oflow.fill(x, y, w);
      } else {
        bins[rx][ry].fill(x, y, w);
      }
    } else if (x < bins[lx][0].xmin || y < bins[0][ly].ymin) {
      uflow.fill(x, y, w);
    } else {
      bins[lx][ly].fill(x, y, w);
    }

    total.fill(x, y, w);
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

  void write(const std::string& filename) const {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);
    file << to_string();
    file.close();
  }
};

#endif  // histogram_h_