#ifndef THRUST_H_
#define THRUST_H_

#include "eshape.h"
#include "histogram.h"

class Thrust : public EShape {
 public:
  Thrust();
  std::tuple<double, std::vector<double>> CalculateThrust(const Event& ev);
};

class TAnalysis {
 public:
  TAnalysis();

  Thrust thr;
  Histo1D thr_hist;
  double wtot;

  void Analyze(const Event& ev);
  void Finalize(const std::string& filename);
};

#endif  // THRUST_H_