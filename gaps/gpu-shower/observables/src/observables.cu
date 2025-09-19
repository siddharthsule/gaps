#include <fstream>

#include "observables.cuh"

// -----------------------------------------------------------------------------
// validate events before binning

__global__ void validate_events(event* events, int* invalid, int n) {
  /**
   * @brief Analyze the event and store the results
   *
   * This function first validates the event, then calculates the observables
   * for the corresponding process. The results are then stored in the
   * histograms and written to a file.
   *
   * @param ev The event object
   * @param invalid The number of invalid events
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------

  event& ev = events[idx];
  ev.set_validity(ev.validate());

  if (!ev.get_validity()) {
    // printf("invalid event\n");
    atomicAdd(invalid, 1);
  }
}

// -----------------------------------------------------------------------------
// analysis

// fill the_histograms (atomically!)
__global__ void fill_histos(analysis* an, const event* events, double* results,
                            int process, int n) {
  /**
   * @brief Fill the histograms with the results of the analysis
   *
   * @param an The analysis object
   * @param ev The event object
   * @param results The array to store the results
   * @param process The process number
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------

  const event& ev = events[idx];

  // LEP: e+ e- -> q qbar
  if (process == 1) {
    // I know forloop can be done here, but just to note what is happening
    an->hists[0].fill(results[20 * idx + 0], ev.get_dxs());  // log10y23
    an->hists[1].fill(results[20 * idx + 1], ev.get_dxs());  // log10y34
    an->hists[2].fill(results[20 * idx + 2], ev.get_dxs());  // log10y45
    an->hists[3].fill(results[20 * idx + 3], ev.get_dxs());  // log10y56
    an->hists[4].fill(results[20 * idx + 4], ev.get_dxs());  // tvalue
    an->hists[5].fill(results[20 * idx + 5], ev.get_dxs());  // tzoomd
    an->hists[6].fill(results[20 * idx + 6], ev.get_dxs());  // hjm
    an->hists[7].fill(results[20 * idx + 7], ev.get_dxs());  // ljm
    an->hists[8].fill(results[20 * idx + 8], ev.get_dxs());  // wjb
    an->hists[9].fill(results[20 * idx + 9], ev.get_dxs());  // njb
  }

  // LHC: p p -> e+ e-
  else if (process == 2) {
    an->hists[0].fill(results[20 * idx + 0], ev.get_dxs());      // zmass
    an->hists[1].fill(results[20 * idx + 1], ev.get_dxs());      // zpt
    an->hists[2].fill(results[20 * idx + 1], ev.get_dxs());      // zptfull
    an->hists[3].fill(results[20 * idx + 2], ev.get_dxs());      // zphi
    an->hists[4].fill(results[20 * idx + 3], ev.get_dxs());      // zrap
    an->hists[5].fill(results[20 * idx + 4], ev.get_dxs() / 2);  // leptpt
    an->hists[5].fill(results[20 * idx + 5], ev.get_dxs() / 2);
    an->hists[6].fill(results[20 * idx + 6], ev.get_dxs() / 2);  // lepteta
    an->hists[6].fill(results[20 * idx + 7], ev.get_dxs() / 2);

    // flavour splitting distribution
    an->hists[7].fill(results[20 * idx + 8], ev.get_dxs());

    // jets
    an->hists[8].fill(results[20 * idx + 9], ev.get_dxs());
    an->hists[9].fill(results[20 * idx + 10], ev.get_dxs());
    an->hists[10].fill(results[20 * idx + 11], ev.get_dxs());
    an->hists[11].fill(results[20 * idx + 12], ev.get_dxs());
    an->hists[12].fill(results[20 * idx + 13], ev.get_dxs());
    an->hists[13].fill(results[20 * idx + 14], ev.get_dxs());
    an->hists[14].fill(results[20 * idx + 15], ev.get_dxs());
    an->hists[15].fill(results[20 * idx + 16], ev.get_dxs());
    an->hists[16].fill(results[20 * idx + 17], ev.get_dxs());
    an->hists[17].fill(results[20 * idx + 18], ev.get_dxs());
    an->hists[18].fill(results[20 * idx + 19], ev.get_dxs());
    an->hists[19].fill(results[20 * idx + 20], ev.get_dxs());

    // Forward-Backward Asymmetry testing
    // an->hists[3].fill(-results[20 * idx + 3], -ev.get_dxs());
  }

  atomicAdd(&an->wtot, ev.get_dxs());
  atomicAdd(&an->ntot, 1.);
}

// -----------------------------------------------------------------------------
// Write Cross Section to file

void write_xsec(double xsec, double xsec_err, const std::string& filename) {
  /**
   * @brief Write the cross-section to a string in YODA format
   *
   * @param xsec The cross-section value
   * @param xsec_err The cross-section error
   * @return std::string The formatted string
   */

  std::stringstream ss;
  ss << "BEGIN YODA_SCATTER1D /_XSEC\n";
  ss << "ErrorBreakdown: {0: {\"\": {up: " << xsec_err << ", dn: " << xsec_err
     << "}}}\n";
  ss << "Path: /_XSEC\n";
  ss << "Title: ~\n";
  ss << "Type: Scatter1D\n";
  ss << "---\n";
  ss << "# xval\t xerr-\t xerr+\t\n";
  ss << std::scientific << std::setprecision(6) << xsec << "\t" << xsec_err
     << "\t" << xsec_err << "\n";
  ss << "END YODA_SCATTER1D_V2\n\n";

  // Write the string to the specified file
  std::ofstream file(filename, std::ios::app);
  if (file.is_open()) {
    file << ss.str();
    file.close();
  }
}

// -----------------------------------------------------------------------------
// run the above kernels
void do_analysis(thrust::device_vector<event>& dv_events, std::string filename,
                 int process, int blocks, int threads) {
  /**
   * @brief Run the analysis
   *
   * @param dv_events device vector of event records
   * @param filename output file name
   * @param process LHC, DIS or LEP
   */

  // device analysis object
  analysis *h_an, *d_an;

  // allocate memory for the device analysis object
  h_an = new analysis(process);
  cudaMalloc(&d_an, sizeof(analysis));
  cudaMemcpy(d_an, h_an, sizeof(analysis), cudaMemcpyHostToDevice);

  // get event data
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // validate the events
  int* d_invalid;
  cudaMalloc(&d_invalid, sizeof(int));
  cudaMemset(d_invalid, 0, sizeof(int));

  validate_events<<<blocks, threads>>>(d_events, d_invalid, n_events);
  sync_gpu_and_check("validate_events");

  int h_invalid;
  cudaMemcpy(&h_invalid, d_invalid, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_invalid);

  if (h_invalid > 0) {
    std::cout << "" << std::endl;
    std::cout << "error: invalid events found" << std::endl;
    std::cout << "number of invalid events: " << h_invalid << "\n";
  }
  std::cout << "" << std::endl;

  // do the analysis

  // Make a vector to store the results of the analysis
  // For now, we set a size 10 per event and use that to store the results
  thrust::device_vector<double> dv_results(20 * n_events, -50.);
  double* d_results = thrust::raw_pointer_cast(dv_results.data());

  if (process == 1) {
    cluster_durham<<<blocks, threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("cluster_durham");

    calculate_ev_shapes<<<blocks, threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("calculate_ev_shapes");
  }

  else if (process == 2) {
    calculate_mczinc<<<blocks, threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("calculate_mczinc");

    cluster_genkt<<<blocks, threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("cluster_genkt");
  }

  // do the analysis
  fill_histos<<<blocks, threads>>>(d_an, d_events, d_results, process,
                                   n_events);
  sync_gpu_and_check("fill_histos");

  // copy the results back to the host
  cudaMemcpy(h_an, d_an, sizeof(analysis), cudaMemcpyDeviceToHost);

  // remove existing file
  std::remove(filename.c_str());

  // Calculate cross-section
  double xsec = h_an->wtot / h_an->ntot;      // Cross-section in pb
  double xsec_err = xsec / sqrt(h_an->ntot);  // Statistical error

  // Print out the total cross-section
  printf("Total cross-section: %.2e nb\n", xsec / 1000.);

  // Write cross-section in YODA format
  write_xsec(xsec, xsec_err, filename);

  for (auto& hist : h_an->hists) {
    if (hist.name[0] != 'h') {
      hist.scale_w(1. / h_an->wtot);
      write(hist, filename);
    }
  }

  // clean up
  delete h_an;
  cudaFree(d_an);
}