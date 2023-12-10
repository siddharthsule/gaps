#ifndef DURHAM_H_
#define DURHAM_H_

#include "histogram.h"
#include "particle.cuh"

// Durham class
class Durham {
 public:
  Durham();
  double Yij(const Vec4& p, const Vec4& q);
  std::vector<double> Cluster(const Event& ev);

  double ecm2;  // square of the center-of-mass energy
};

// DAnalysis class
class DAnalysis {
 public:
  DAnalysis();

  Durham dur;
  std::vector<Histo1D> dur_hists;
  double wtot;

  void Analyze(const Event& ev);
  void Finalize(const std::string& filename);
};

#endif  // DURHAM_H_