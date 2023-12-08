#include "histogram.h"

#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>

Bin1D::Bin1D(double xmin, double xmax)
    : xmin(xmin), xmax(xmax), w(0.0), w2(0.0), wx(0.0), wx2(0.0), n(0.0) {}

std::string Bin1D::Format(const std::string& tag) const {
  std::stringstream ss;
  ss << std::scientific << std::setprecision(6);
  ss << tag << "\t" << tag << "\t" << w << "\t" << w2 << "\t" << wx << "\t" 
     << wx2 << "\t" << static_cast<int>(n);
  return ss.str();
}

std::string Bin1D::ToString() const {
  std::stringstream ss;
  ss << std::scientific << std::setprecision(6);
  ss << xmin << "\t" << xmax << "\t" << w << "\t" << w2 << "\t" << wx << "\t" 
     << wx2 << "\t" << static_cast<int>(n);
  return ss.str();
}

double Bin1D::Width() const {
  return xmax - xmin;
}

void Bin1D::Fill(double x, double weight) {
  this->w += weight;
  w2 += weight * weight;
  wx += weight * x;
  wx2 += weight * weight * x;
  n += 1.0;
}

void Bin1D::ScaleW(double scale) {
  w *= scale;
  w2 *= scale * scale;
  wx *= scale;
  wx2 *= scale * scale;
}


// Constructor for Histo1D
Histo1D::Histo1D(int nbin, double xmin, double xmax, const std::string& name)
    : name(name), uflow(xmin - 100., xmin), oflow(xmax, xmax + 100.),
      total(xmin - 100., xmax + 100.), scale(1.) {
  double width = (xmax - xmin) / nbin;
  for (int i = 0; i < nbin; ++i) {
    double xlow = xmin + i * width;
    double xhigh = xlow + width;
    bins.push_back(Bin1D(xlow, xhigh));
  }
}

std::string Histo1D::ToString() const {
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

void Histo1D::Fill(double x, double w) {
  if (x < uflow.xmax) {
    uflow.Fill(x, w);
  } else if (x >= oflow.xmin) {
    oflow.Fill(x, w);
  } else {
    int bin = static_cast<int>((x - bins[0].xmin) / bins[0].Width());
    bins[bin].Fill(x, w);
  }
  total.Fill(x, w);
}

void Histo1D::ScaleW(double scale) {
  for (auto& bin : bins) {
    bin.ScaleW(scale);
  }
  uflow.ScaleW(scale);
  oflow.ScaleW(scale);
  total.ScaleW(scale);
  this->scale *= scale;
}

void Histo1D::Write(const std::string& filename) const {
  std::ofstream file;
  file.open(filename, std::ios::out | std::ios::app);
  file << ToString();
  file.close();
}

