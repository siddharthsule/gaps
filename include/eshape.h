#ifndef ESHAPE_H_
#define ESHAPE_H_

#include <vector>

#include "particle.cuh"
#include "histogram.h"

// EShape class
class EShape {
 public:
  EShape();
  std::vector<Vec4> GetMomenta(const Event& ev);

  // Thrust
  std::tuple<double, Vec4> CalculateThrust(const Event& ev);

  // Jet Mass and Broadening
  std::vector<double> CalculateJetMBr(const Event& ev);

};

// EAnalysis class
class EAnalysis {
 public:
  EAnalysis();

  EShape eshape;

  // 0 = Thrust, 1 = Thrust, Zoomed In
  // 2 = Heavy Jet Mass, 3 = Light Jet Mass
  // 4 = Wide Jet Broadening, 5 = Narrow Jet Broadening
  std::vector<Histo1D> eshape_hists;
  double wtot;

  void Analyze(const Event& ev);
  void Finalize(const std::string& filename);
};

#endif  // ESHAPE_H_