#include <fstream>

#include "observables.cuh"

// -----------------------------------------------------------------------------
// Validate Events before binning

__global__ void validateEvents(Event* events, int* invalid, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];
  ev.SetValidity(ev.Validate());

  if (!ev.GetValidity()) {
    // printf("Invalid Event\n");
    atomicAdd(invalid, 1);
  }
}

// -----------------------------------------------------------------------------
// Analysis

// Fill theHistograms (Atomically!)
__global__ void fillHistos(Analysis* an, Event* events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];

  an->hists[0].Fill(ev.GetY23(), ev.GetDxs());
  an->hists[1].Fill(ev.GetY34(), ev.GetDxs());
  an->hists[2].Fill(ev.GetY45(), ev.GetDxs());
  an->hists[3].Fill(ev.GetY56(), ev.GetDxs());
  an->hists[4].Fill(ev.GetThr(), ev.GetDxs());
  an->hists[5].Fill(ev.GetThr(), ev.GetDxs());
  an->hists[6].Fill(ev.GetHJM(), ev.GetDxs());
  an->hists[7].Fill(ev.GetLJM(), ev.GetDxs());
  an->hists[8].Fill(ev.GetWJB(), ev.GetDxs());
  an->hists[9].Fill(ev.GetNJB(), ev.GetDxs());

  // Dalitz Plot is OFF
  // an->dalitz.Fill(ev.GetDalitz(0), ev.GetDalitz(1), ev.GetDxs());

  atomicAdd(&an->wtot, ev.GetDxs());
  atomicAdd(&an->ntot, 1.);
}

// Run the above kernels
void doAnalysis(thrust::device_vector<Event>& d_events, std::string filename) {
  /**
   * Only Place for a Host Object
   * ----------------------------
   *
   * While most of the work is done on the device, one cannot directly write
   * to a file from the device. Therefore, we will create a host object to
   * store the histograms, and then copy the results back to the host for
   * writing to file.
   */

  // Device Analysis Object
  Analysis *h_an, *d_an;

  // Allocate memory for the device analysis object
  h_an = new Analysis();
  cudaMalloc(&d_an, sizeof(Analysis));
  cudaMemcpy(d_an, h_an, sizeof(Analysis), cudaMemcpyHostToDevice);

  // Get Event Data
  int N = d_events.size();
  Event* d_events_ptr = thrust::raw_pointer_cast(d_events.data());

  // Validate the Events
  int* d_invalid;
  cudaMalloc(&d_invalid, sizeof(int));
  cudaMemset(d_invalid, 0, sizeof(int));

  validateEvents<<<(N + 255) / 256, 256>>>(d_events_ptr, d_invalid, N);
  syncGPUAndCheck("validateEvents");

  int h_invalid;
  cudaMemcpy(&h_invalid, d_invalid, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_invalid);

  if (h_invalid > 0) {
    std::cout << "" << std::endl;
    std::cout << "ERROR: Invalid Events Found" << std::endl;
    std::cout << "Number of Invalid Events: " << h_invalid << "\n";
  }

  // Calculare the Observables
  doCluster<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("doCluster");

  calculateThr<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("calculateThr");

  calculateJetMBr<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("calculateJetMBr");

  /**
   * Why is the Dalitz Plot off?
   * ---------------------------
   *
   * While the Dalitz analysis also benefits from the GPU parallelisation, the
   * writing of the data to file severely limits the performance, as instead of
   * the usual 100 bins, we have 100^2 = 1000 bins. This takes around 0.04s,
   * which is minute in the C++ case, but is in fact 40% of the total analysis
   * time! So for our tests, we keep this off
   *
   * If you want to turn it on, uncomment the lines in this file, and it's
   * equivalent in the 'observables.cpp' file.
   */
  // calculateDalitz<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  // syncGPUAndCheck("calculateDalitz");

  // Do the Analysis
  fillHistos<<<(N + 255) / 256, 256>>>(d_an, d_events_ptr, N);
  syncGPUAndCheck("fillHistos");

  // Copy the results back to the host
  cudaMemcpy(h_an, d_an, sizeof(Analysis), cudaMemcpyDeviceToHost);

  // Normalize the histograms
  for (auto& hist : h_an->hists) {
    hist.ScaleW(1. / h_an->ntot);
  }

  // Dalitz Plot is OFF
  // h_an->dalitz.ScaleW(1. / h_an->ntot);

  // Remove existing file
  std::remove(filename.c_str());

  // Write the histograms to file
  Write(h_an->hists[0], "/gaps/log10y23\n", filename);
  Write(h_an->hists[1], "/gaps/log10y34\n", filename);
  Write(h_an->hists[2], "/gaps/log10y45\n", filename);
  Write(h_an->hists[3], "/gaps/log10y56\n", filename);
  Write(h_an->hists[4], "/gaps/tvalue\n", filename);
  Write(h_an->hists[5], "/gaps/tzoomd\n", filename);
  Write(h_an->hists[6], "/gaps/hjm\n", filename);
  Write(h_an->hists[7], "/gaps/ljm\n", filename);
  Write(h_an->hists[8], "/gaps/wjb\n", filename);
  Write(h_an->hists[9], "/gaps/njb\n", filename);

  // Dalitz Plot is OFF
  // Write(h_an->dalitz, "/gaps/dalitz\n", filename);

  // Clean up
  delete h_an;
  cudaFree(d_an);
}