#ifndef durham_h_
#define durham_h_

#include "dalitz.h"
#include "event.h"
#include "eventshapes.h"
#include "histogram.h"
#include "jetrates.h"

/**
 * slight difference between c++ and cuda codes
 * --------------------------------------------
 *
 * here, we don't need analyse and finalize functions to be outside the class.
 * so we can just include them in the class definition. also we have to use the
 * same class for all events!
 */
class analysis {
 public:
  histo1d hists[10];
  histo2d dalitz;

  double wtot;  // scale by weight for 1/sigma d(sigma)/d observable
  double ntot;  // scale by number for d(sigma)/d observable

 public:
  analysis() : wtot(0.), ntot(0.) {
    hists[0] = histo1d(-4.3, -0.3, "/gaps/log10y23\n");
    hists[1] = histo1d(-4.3, -0.3, "/gaps/log10y34\n");
    hists[2] = histo1d(-4.3, -0.3, "/gaps/log10y45\n");
    hists[3] = histo1d(-4.3, -0.3, "/gaps/log10y56\n");
    hists[4] = histo1d(0., 0.5, "/gaps/tvalue\n");
    hists[5] = histo1d(0., 0.5, "/gaps/tzoomd\n");
    hists[6] = histo1d(0., 1., "/gaps/hjm\n");
    hists[7] = histo1d(0., 0.5, "/gaps/ljm\n");
    hists[8] = histo1d(0., 0.5, "/gaps/wjb\n");
    hists[9] = histo1d(0., 0.2, "/gaps/njb\n");

    dalitz = histo2d(0., 1., 0., 1., "/gaps/dalitz\n");
  }

  // member functions here, but not in cuda
  void analyze(event& ev);
  void finalize(const std::string& filename);
};

#endif  // durham_h_