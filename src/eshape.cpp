#include "eshape.h"

#include <algorithm>
#include <cmath>
#include <vector>

EShape::EShape() {}

std::vector<std::vector<double>> EShape::GetMomenta(const Event& ev) {
  std::vector<std::vector<double>> moms;

  // Get Momenta for particles starting at 2
  for (size_t i = 2; i < ev.GetSize(); ++i) {
    Vec4 pmom = ev.GetParton(i).GetMom();
    std::vector<double> mom = {pmom[1], pmom[2], pmom[3]};
    moms.push_back(mom);
  }

  // Sort by mag2 in descending order
  std::sort(moms.begin(), moms.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
    return MagnitudeSquared(a) > MagnitudeSquared(b);
  });

  return moms;
}

double Magnitude(const std::vector<double>& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

double MagnitudeSquared(const std::vector<double>& v) {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

std::vector<double> Add(const std::vector<double>& a, const std::vector<double>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

std::vector<double> Subtract(const std::vector<double>& a, const std::vector<double>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::vector<double> Multiply(const std::vector<double>& v, const double& s) {
    return {v[0] * s, v[1] * s, v[2] * s};
}

double DotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::vector<double> CrossProduct(const std::vector<double>& a, const std::vector<double>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}