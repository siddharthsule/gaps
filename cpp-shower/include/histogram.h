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
      : xmin(xmin), xmax(xmax), w(0.0), w2(0.0), wx(0.0), wx2(0.0), n(0.0) {}

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
    n += 1.0;
  }

  void ScaleW(double scale) {
    w *= scale;
    w2 *= scale * scale;
    wx *= scale;
    wx2 *= scale * scale;
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
  Histo1D(double xmin = 0.0, double xmax = 1.0, const std::string& name = "hst")
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

#endif  // HISTOGRAM_H_