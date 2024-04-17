#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "base.h"

// Bin1D class
class Bin1D {
 public:
  double xmin, xmax, w, w2, wx, wx2, n;

 public:
  Bin1D(double xmin, double xmax)
      : xmin(xmin), xmax(xmax), w(0.), w2(0.), wx(0.), wx2(0.), n(0.) {}

  std::string Format(const std::string& tag) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << tag << "\t" << tag << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << xmin << "\t" << xmax << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  double Width() const { return xmax - xmin; }

  void Fill(double x, double weight) {
    this->w += weight;
    w2 += weight * weight;
    wx += weight * x;
    wx2 += weight * weight * x;
    n += 1.;
  }

  void ScaleW(double scale) {
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
  Bin2D(double xmin, double xmax, double ymin, double ymax)
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

  std::string Format(const std::string& tag) const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << tag << "\t" << tag << "\t" << w << "\t" << w2 << "\t" << wx << "\t"
       << wx2 << "\t" << wy << "\t" << wy2 << "\t" << static_cast<int>(n);
    return ss.str();
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(6);
    ss << xmin << "\t" << xmax << "\t" << ymin << "\t" << ymax << "\t" << w
       << "\t" << w2 << "\t" << wx << "\t" << wx2 << "\t" << wy << "\t" << wy2
       << "\t" << wxy << "\t" << static_cast<int>(n);
    return ss.str();
  }

  double WidthX() const { return xmax - xmin; }
  double WidthY() const { return ymax - ymin; }

  void Fill(double x, double y, double weight) {
    this->w += weight;
    w2 += weight * weight;
    wx += weight * x;
    wx2 += weight * weight * x;
    wy += weight * y;
    wy2 += weight * weight * y;
    wxy += weight * x * y;
    n += 1.;
  }

  void ScaleW(double scale) {
    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
    wy *= scale;
    wy2 *= scale * scale;
    wxy *= scale * scale;
  }
};

// Histo1D class
class Histo1D {
 public:
  std::string name;
  std::vector<Bin1D> bins;
  Bin1D uflow;
  Bin1D oflow;
  Bin1D total;
  double scale;

 public:
  // Constructor for Histo1D
  Histo1D(double xmin = 0., double xmax = 1., const std::string& name = "hst")
      : name(name),
        uflow(xmin - 100., xmin),
        oflow(xmax, xmax + 100.),
        total(xmin - 100., xmax + 100.),
        scale(1.) {
    double width = (xmax - xmin) / nBins;
    for (int i = 0; i < nBins; ++i) {
      double xlow = xmin + i * width;
      double xhigh = xlow + width;
      bins.push_back(Bin1D(xlow, xhigh));
    }
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << "BEGIN YODA_HISTO1D " << name << "\n\n";
    ss << "Path=" << name << "\n\n";
    ss << "ScaledBy=" << scale << "\n";
    ss << "Title=\nType=Histo1D\n";
    ss << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    ss << total.Format("Total") << "\n";
    ss << uflow.Format("Underflow") << "\n";
    ss << oflow.Format("Overflow") << "\n";
    ss << "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
    for (const auto& bin : bins) {
      ss << bin.ToString() << "\n";
    }
    ss << "END YODA_HISTO1D\n\n";
    return ss.str();
  }

  void Fill(double x, double w) {
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
        oflow.Fill(x, w);
      } else {
        bins[r].Fill(x, w);
      }
    } else if (x < bins[l].xmin) {
      uflow.Fill(x, w);
    } else {
      bins[l].Fill(x, w);
    }

    total.Fill(x, w);
  }

  void ScaleW(double scale) {
    for (auto& bin : bins) {
      bin.ScaleW(scale);
    }
    uflow.ScaleW(scale);
    oflow.ScaleW(scale);
    total.ScaleW(scale);
    this->scale *= scale;
  }

  void Write(const std::string& filename) const {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);
    file << ToString();
    file.close();
  }
};

// Histo2D class
class Histo2D {
public:
  std::string name;
  std::vector<std::vector<Bin2D>> bins;
  Bin2D uflow;
  Bin2D oflow;
  Bin2D total;
  double scale;

public:
  // Constructor for Histo2D
  Histo2D(double xmin = 0., double xmax = 1., double ymin = 0., double ymax = 1., const std::string& name = "hst")
    : name(name),
      uflow(xmin - 100., xmin, ymin - 100., ymin),
      oflow(xmax, xmax + 100., ymax, ymax + 100.),
      total(xmin - 100., xmax + 100., ymin - 100., ymax + 100.),
      scale(1.) {
    double xwidth = (xmax - xmin) / nBins;
    double ywidth = (ymax - ymin) / nBins;
    for (int i = 0; i < nBins; ++i) {
      double xlow = xmin + i * xwidth;
      double xhigh = xlow + xwidth;
      std::vector<Bin2D> binRow;
      for (int j = 0; j < nBins; ++j) {
        double ylow = ymin + j * ywidth;
        double yhigh = ylow + ywidth;
        binRow.push_back(Bin2D(xlow, xhigh, ylow, yhigh));
      }
      bins.push_back(binRow);
    }
  }

  std::string ToString() const {
    std::stringstream ss;
    ss << "BEGIN YODA_HISTO2D " << name << "\n\n";
    ss << "Path=" << name << "\n\n";
    ss << "ScaledBy=" << scale << "\n";
    ss << "Title=\nType=Histo2D\n";
    ss << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2\tnumEntries\n";
    ss << total.Format("Total") << "\n";
    //ss << uflow.Format("Underflow") << "\n";
    //ss << oflow.Format("Overflow") << "\n";
    ss << "# xlow\txhigh\tylow\tyhigh\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2\tnumEntries\n";
    for (const auto& binRow : bins) {
      for (const auto& bin : binRow) {
        ss << bin.ToString() << "\n";
      }
    }
    ss << "END YODA_HISTO2D\n\n";
    return ss.str();
  }

  void Fill(double x, double y, double w) {
    // Find the bin for the x-coordinate
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

    // Find the bin for the y-coordinate
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

    // Fill the appropriate bin
    if (x > bins[rx][0].xmin && y > bins[0][ry].ymin) {
      if (x > bins[rx][0].xmax || y > bins[0][ry].ymax) {
        oflow.Fill(x, y, w);
      } else {
        bins[rx][ry].Fill(x, y, w);
      }
    } else if (x < bins[lx][0].xmin || y < bins[0][ly].ymin) {
      uflow.Fill(x, y, w);
    } else {
      bins[lx][ly].Fill(x, y, w);
    }

    total.Fill(x, y, w);
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

  void Write(const std::string& filename) const {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);
    file << ToString();
    file.close();
  }
};

#endif  // HISTOGRAM_H_