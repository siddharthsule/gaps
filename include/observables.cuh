#ifndef DURHAM_CUH_
#define DURHAM_CUH_

// Histogram and Parton include relevant headers
#include "histogram.cuh"
#include "parton.cuh"

/**
 * Observables and their Anaylsis
 * -------------------------------
 *
 * This file contains the necessary classes and functions to perform a Durham
 * algorithm analysis on a set of events, as well as thrust and Jet massses +
 * broadenings. The observables are calculated here and then analysed using
 * Rivet[1]
 *
 * [1] Rivet: https://rivet.hepforge.org/
 */

// Jet rates using Durham algorithm
__device__ double Yij(const Vec4& p, const Vec4& q, double ecm2);
__global__ void doCluster(Event* events, int N);

// Event Shapes
__device__ void bubbleSort(Vec4* moms, int n);
__global__ void calculateThr(Event* events, int N);
__global__ void calculateJetMBr(Event* events, int N);

// Dalitz Plot
__global__ void calculateDalitz(Event* events, int N);

class Analysis {
 public:
  // Similar to Histo1D in C++/Rivet, just split into Host / Device Components
  Histo1D hists[10];
  Histo2D dalitz;

  double wtot;  // Scale by Weight for 1/sigma d(sigma)/d Observable
  double ntot;  // Scale by Number for d(sigma)/d Observable

  // Can't have strings as device variables, in future could use char*
  __host__ __device__ Analysis() : wtot(0.0), ntot(0.0) {
    hists[0] = Histo1D(-4.3, -0.3);  // /gaps/log10y23
    hists[1] = Histo1D(-4.3, -0.3);  // /gaps/log10y34
    hists[2] = Histo1D(-4.3, -0.3);  // /gaps/log10y45
    hists[3] = Histo1D(-4.3, -0.3);  // /gaps/log10y56
    hists[4] = Histo1D(0.0, 0.5);    // "/gaps/tvalue"
    hists[5] = Histo1D(0.0, 0.5);    // "/gaps/tzoomd"
    hists[6] = Histo1D(0.0, 1.0);    // "/gaps/hjm"
    hists[7] = Histo1D(0.0, 0.5);    // "/gaps/ljm"
    hists[8] = Histo1D(0.0, 0.5);    // "/gaps/wjb"
    hists[9] = Histo1D(0.0, 0.2);    // "/gaps/njb"

    dalitz = Histo2D(0.0, 1.0, 0.0, 1.0);
  }
};

// Validate Events - Colour and Momentum Conservation
__global__ void validateEvents(Event* events, int N);

// Fill Histograms simultaneously
__global__ void fillHistos(Analysis* an, Event* events, int N);

// Analysis wrapped in a function
void doAnalysis(thrust::device_vector<Event>& d_events, std::string filename);

#endif  // DURHAM_CUH_