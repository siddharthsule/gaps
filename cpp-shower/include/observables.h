#ifndef DURHAM_H_
#define DURHAM_H_

#include "histogram.h"
#include "parton.h"

#include <algorithm>

// Jet Rates using the Durham Algorithm
double Yij(const Vec4& p, const Vec4& q);
void Cluster(Event& ev);

// Event Shapes
std::vector<Vec4> GetMomenta(const Event& ev);
void CalculateThrust(Event& ev);
void CalculateJetMBr(Event& ev);

/**
 * Slight Difference between C++ and CUDA codes
 * --------------------------------------------
 *
 * Here, we don't need Analyse and Finalize functions to be outside the class.
 * So we can just include them in the class definition. Also we have to use the
 * same class for all events!
 */
class Analysis {
  std::vector<Histo1D> hists;
  double wtot;
  double ntot;

 public:
  Analysis() : wtot(0.0), ntot(0.0) {
    hists.push_back(Histo1D(-4.3, -0.3, "/gaps/log10y23\n"));
    hists.push_back(Histo1D(-4.3, -0.3, "/gaps/log10y34\n"));
    hists.push_back(Histo1D(-4.3, -0.3, "/gaps/log10y45\n"));
    hists.push_back(Histo1D(-4.3, -0.3, "/gaps/log10y56\n"));
    hists.push_back(Histo1D(0.0, 0.5, "/gaps/tvalue\n"));
    hists.push_back(Histo1D(0.0, 0.5, "/gaps/tzoomd\n"));
    hists.push_back(Histo1D(0.0, 1.0, "/gaps/hjm\n"));
    hists.push_back(Histo1D(0.0, 0.5, "/gaps/ljm\n"));
    hists.push_back(Histo1D(0.0, 0.5, "/gaps/wjb\n"));
    hists.push_back(Histo1D(0.0, 0.2, "/gaps/njb\n"));
  }

  void Analyze(Event& ev);
  void Finalize(const std::string& filename);
};

#endif  // DURHAM_H_