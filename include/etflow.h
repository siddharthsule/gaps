#ifndef ETFLOW_H_
#define ETFLOW_H_

#include "eshape.h"
#include "thrust.h"
#include "histogram.h"

// ETFlow class
class ETFlow : public EShape {
 public:
  ETFlow();
  std::vector<double> CalculateETFlow(const Event& ev);

  Thrust thr;
};

// ETAnalysis class
class ETAnalysis {
 public:
  ETAnalysis();

  ETFlow etflow;
  std::vector<Histo1D> etflow_hists;
  double wtot;

  void Analyze(const Event& ev);
  void Finalize(const std::string& filename);

};

#endif  // ETFLOW_H_