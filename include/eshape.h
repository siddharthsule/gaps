#ifndef ESHAPE_H_
#define ESHAPE_H_

#include <vector>
#include "particle.cuh"

// EShape class
class EShape {
 public:
  EShape();
  std::vector<std::vector<double>> GetMomenta(const Event& ev);

};

double Magnitude(const std::vector<double>& v);
double MagnitudeSquared(const std::vector<double>& v);
std::vector<double> Add(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> Subtract(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> Multiply(const std::vector<double>& v, const double& s);
double DotProduct(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> CrossProduct(const std::vector<double>& a, const std::vector<double>& b);

#endif  // ESHAPE_H_