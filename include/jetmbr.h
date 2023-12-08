#ifndef JETMBR_H_
#define JETMBR_H_

#include "eshape.h"
#include "thrust.h"
#include "histogram.h"

// JetMBr class
class JetMBr : public EShape {
 public:
  JetMBr();
  std::vector<double> CalculateJetMBr(const Event& ev);

  Thrust thr;
};

// MBrAnalysis class
class MBrAnalysis {
 public:
  MBrAnalysis();

  JetMBr jetmbr;
  std::vector<Histo1D> jetmbr_hists;
  double wtot;

  void Analyze(const Event& ev);
  void Finalize(const std::string& filename);
};

#endif  // JETMBR_H_