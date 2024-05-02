#ifndef durham_cuh_
#define durham_cuh_

// histogram and parton include relevant headers
#include "event.cuh"
#include "histogram.cuh"

// add analyses here
#include "dalitz.cuh"
#include "eventshapes.cuh"
#include "jetrates.cuh"

/**
 * observables and their anaylsis
 * -------------------------------
 *
 * this file contains the necessary classes and functions to perform a durham
 * algorithm analysis on a set of events, as well as thrust and jet massses +
 * broadenings. the observables are calculated here and then analysed using
 * rivet[1]
 *
 * [1] rivet: https://rivet.hepforge.org/
 */

class analysis {
 public:
  // similar to histo1d in c++/rivet, just split into host / device components
  histo1d hists[10];
  histo2d dalitz;

  double wtot;  // scale by weight for 1/sigma d(sigma)/d observable
  double ntot;  // scale by number for d(sigma)/d observable

  // can't have strings as device variables, in future could use char*
  __host__ __device__ analysis() : wtot(0.), ntot(0.) {
    hists[0] = histo1d(-4.3, -0.3);  // /gaps/log10y23
    hists[1] = histo1d(-4.3, -0.3);  // /gaps/log10y34
    hists[2] = histo1d(-4.3, -0.3);  // /gaps/log10y45
    hists[3] = histo1d(-4.3, -0.3);  // /gaps/log10y56
    hists[4] = histo1d(0., 0.5);     // "/gaps/tvalue"
    hists[5] = histo1d(0., 0.5);     // "/gaps/tzoomd"
    hists[6] = histo1d(0., 1.);      // "/gaps/hjm"
    hists[7] = histo1d(0., 0.5);     // "/gaps/ljm"
    hists[8] = histo1d(0., 0.5);     // "/gaps/wjb"
    hists[9] = histo1d(0., 0.2);     // "/gaps/njb"

    dalitz = histo2d(0., 1., 0., 1.);  // '/gaps/dalitz"
  }
};

// validate events - colour and momentum conservation
__global__ void validate_events(event* events, int* invalid, int n);

// fill histograms simultaneously
__global__ void fill_histos(analysis* an, event* events, int n);

// analysis wrapped in a function
void do_analysis(thrust::device_vector<event>& d_events, std::string filename);

#endif  // durham_cuh_