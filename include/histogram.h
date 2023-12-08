#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <string>
#include <vector>

// Bin1D class
class Bin1D {
 public:
  Bin1D(double xmin, double xmax);

  std::string Format(const std::string& tag) const;
  std::string ToString() const;
  
  double Width() const;
  void Fill(double x, double w);
  void ScaleW(double scale);

  double xmin;
  double xmax;
  double w;
  double w2;
  double wx;
  double wx2;
  double n;
};

// Histo1D class
class Histo1D {
 public:
  Histo1D(int nbin, double xmin, double xmax, const std::string& name);

  std::string ToString() const;
  void Fill(double x, double w);
  void ScaleW(double scale);
  void Write(const std::string& filename) const;

  std::string name;
  std::vector<Bin1D> bins;
  Bin1D uflow;
  Bin1D oflow;
  Bin1D total;
  double scale;
};

#endif  // HISTOGRAM_H_