#ifndef DURHAM_H_
#define DURHAM_H_

#include "histogram.h"
#include "event.h"

#include "jetrates.h"
#include "eventshapes.h"
#include "dalitz.h"

/**
 * Slight Difference between C++ and CUDA codes
 * --------------------------------------------
 *
 * Here, we don't need Analyse and Finalize functions to be outside the class.
 * So we can just include them in the class definition. Also we have to use the
 * same class for all events!
 */
class Analysis {
 public:
  Histo1D hists[10];
  Histo2D dalitz;

  double wtot;  // Scale by Weight for 1/sigma d(sigma)/d Observable
  double ntot;  // Scale by Number for d(sigma)/d Observable

 public:
  Analysis() : wtot(0.0), ntot(0.0) {
    hists[0] = Histo1D(-4.3, -0.3, "/gaps/log10y23\n");
    hists[1] = Histo1D(-4.3, -0.3, "/gaps/log10y34\n");
    hists[2] = Histo1D(-4.3, -0.3, "/gaps/log10y45\n");
    hists[3] = Histo1D(-4.3, -0.3, "/gaps/log10y56\n");
    hists[4] = Histo1D(0.0, 0.5, "/gaps/tvalue\n");
    hists[5] = Histo1D(0.0, 0.5, "/gaps/tzoomd\n");
    hists[6] = Histo1D(0.0, 1.0, "/gaps/hjm\n");
    hists[7] = Histo1D(0.0, 0.5, "/gaps/ljm\n");
    hists[8] = Histo1D(0.0, 0.5, "/gaps/wjb\n");
    hists[9] = Histo1D(0.0, 0.2, "/gaps/njb\n");

    dalitz = Histo2D(0.0, 1.0, 0.0, 1.0, "/gaps/dalitz\n");
  }

  // Member Functions Here, but not in CUDA
  void Analyze(Event& ev);
  void Finalize(const std::string& filename);
};

#endif  // DURHAM_H_