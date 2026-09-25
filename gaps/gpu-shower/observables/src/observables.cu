#include <fstream>
#include <iomanip>
#include <sstream>

#include "observables.cuh"

// -----------------------------------------------------------------------------
// validate events before binning

__global__ void validate_events(event* events, int* invalid, int* overflowed,
                                int n) {
  /**
   * @brief Validate the event
   *
   * This function checks if the event is valid and sets the validity flag
   * accordingly.
   *
   * @param events The array of event objects
   * @param invalid The number of invalid events
   * @param overflowed How many of those ran out of room in the record
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  event& ev = events[idx];
  // ---------------------------------------------
  ev.set_validity(ev.validate());

  if (!ev.get_validity()) {
    atomicAdd(invalid, 1);
    if (ev.get_overflowed()) atomicAdd(overflowed, 1);
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
   * @param events The events array
   * @param results The array to store the results
   * @param process The process number
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Event Preamble
  const event& ev = events[idx];
  // ---------------------------------------------

  // LEP: e+ e- -> q qbar
  if (process == 1) {
    // fill histograms
    // I know forloop can be done here, but just to note what is happening
    an->hists[0].fill(results[20 * idx + 0], ev.get_dxs());       // log10y23
    an->hists[1].fill(results[20 * idx + 1], ev.get_dxs());       // log10y34
    an->hists[2].fill(results[20 * idx + 2], ev.get_dxs());       // log10y45
    an->hists[3].fill(results[20 * idx + 3], ev.get_dxs());       // log10y56
    an->hists[4].fill(1. - results[20 * idx + 4], ev.get_dxs());  // tvalue
    an->hists[5].fill(1. - results[20 * idx + 5], ev.get_dxs());  // tzoomd
    an->hists[6].fill(results[20 * idx + 6], ev.get_dxs());       // hjm
    an->hists[7].fill(results[20 * idx + 7], ev.get_dxs());       // ljm
    an->hists[8].fill(results[20 * idx + 8], ev.get_dxs());       // wjb
    an->hists[9].fill(results[20 * idx + 9], ev.get_dxs());       // njb
    an->hists[10].fill(ev.get_size() - 2, ev.get_dxs());          // nump
    an->hists[22].fill(results[20 * idx + 11], ev.get_dxs());     // tjb

    // ALEPH
    an->hists[11].fill(results[20 * idx + 4], ev.get_dxs());
    an->hists[12].fill(results[20 * idx + 6] * results[20 * idx + 6],
                       ev.get_dxs());
    an->hists[13].fill(results[20 * idx + 8], ev.get_dxs());
    an->hists[14].fill(results[20 * idx + 12], ev.get_dxs());
    an->hists[15].fill(results[20 * idx + 11], ev.get_dxs());
    an->hists[16].fill(-log(pow(10., results[20 * idx + 0])), ev.get_dxs());
    an->hists[17].fill(-log(pow(10., results[20 * idx + 1])), ev.get_dxs());
    an->hists[18].fill(-log(pow(10., results[20 * idx + 2])), ev.get_dxs());
    an->hists[19].fill(-log(pow(10., results[20 * idx + 3])), ev.get_dxs());

    // L3
    an->hists[20].fill(results[20 * idx + 10], ev.get_dxs());  // L3 nch
    fill_log_scaled_mom(ev, an->hists[21], ev.get_dxs());      // L3 xi
  }

  // LHC: p p -> e+ e-
  else if (process == 2) {
    // fill histograms
    an->hists[0].fill(results[20 * idx + 0], ev.get_dxs());      // zmass
    an->hists[1].fill(results[20 * idx + 1], ev.get_dxs());      // zpt
    an->hists[2].fill(results[20 * idx + 1], ev.get_dxs());      // zptfull
    an->hists[3].fill(results[20 * idx + 2], ev.get_dxs());      // zphi
    an->hists[4].fill(results[20 * idx + 3], ev.get_dxs());      // zrap
    an->hists[5].fill(results[20 * idx + 4], ev.get_dxs() / 2);  // leptpt
    an->hists[5].fill(results[20 * idx + 5], ev.get_dxs() / 2);
    an->hists[6].fill(results[20 * idx + 6], ev.get_dxs() / 2);  // lepteta
    an->hists[6].fill(results[20 * idx + 7], ev.get_dxs() / 2);

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

    // Forward-Backward Asymmetry testing
    // an->hists[3].fill(-results[20 * idx + 3], -ev.get_dxs());
  }

  // weighted total
  atomicAdd(&an->wtot, ev.get_dxs());
  atomicAdd(&an->ntot, 1.);
}

// -----------------------------------------------------------------------------
// Write Cross Section to file

void write_xsec(double xsec, double xsec_err, const std::string& filename) {
  /**
   * @brief Append the cross-section to a file in YODA format
   *
   * @param xsec The cross-section value
   * @param xsec_err The cross-section error
   * @param filename The file to append to
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
void do_analysis(thrust::device_vector<event>& dv_events, const params& p,
                 int blocks) {
  /**
   * @brief Run the analysis
   *
   * @param dv_events device vector of event records
   * @param p run parameters (process, storage_file, threads, ...)
   * @param blocks number of CUDA blocks
   */

  // device analysis object
  analysis *h_an, *d_an;

  // allocate memory for the device analysis object
  h_an = new analysis(p.process);
  cudaMalloc(&d_an, sizeof(analysis));
  cudaMemcpy(d_an, h_an, sizeof(analysis), cudaMemcpyHostToDevice);

  // get event data
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // validate the events

  // Events failing momentum or colour conservation, or overflowed
  int* d_invalid;
  cudaMalloc(&d_invalid, sizeof(int));
  cudaMemset(d_invalid, 0, sizeof(int));

  // Of those, the ones that ran out of room in the record
  int* d_overflowed;
  cudaMalloc(&d_overflowed, sizeof(int));
  cudaMemset(d_overflowed, 0, sizeof(int));

  validate_events<<<blocks, p.threads>>>(d_events, d_invalid, d_overflowed,
                                         n_events);
  sync_gpu_and_check("validate_events");

  int h_invalid, h_overflowed;
  cudaMemcpy(&h_invalid, d_invalid, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(&h_overflowed, d_overflowed, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_invalid);
  cudaFree(d_overflowed);

  if (h_invalid > 0) {
    std::cout << "" << std::endl;
    std::cout << "error: invalid events found" << std::endl;
    std::cout << "number of invalid events: " << h_invalid << "\n";

    if (h_overflowed > 0) {
      std::cout << h_overflowed << " of them ran out of room in the record"
                << std::endl;
      std::cout << "Consider increasing max_particles, default: "
                << max_particles << std::endl;
    }
  }
  std::cout << "" << std::endl;

  // do the analysis

  // Store the results of the analysis, 20 slots per event
  thrust::device_vector<double> dv_results(20 * n_events, -50.);
  double* d_results = thrust::raw_pointer_cast(dv_results.data());

  if (p.process == 1) {
    // jet rates
    cluster_durham<<<blocks, p.threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("cluster_durham");

    // event shapes
    calculate_ev_shapes<<<blocks, p.threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("calculate_ev_shapes");

    // L3
    calculate_chargedmult<<<blocks, p.threads>>>(d_events, d_results, 10,
                                                 n_events);
    sync_gpu_and_check("calculate_chargedmult");
  }

  else if (p.process == 2) {
    // calculate Z boson observables
    calculate_mczinc<<<blocks, p.threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("calculate_mczinc");

    // Gen kt Algorithm
    cluster_genkt<<<blocks, p.threads>>>(d_events, d_results, n_events);
    sync_gpu_and_check("cluster_genkt");
  }

  // do the analysis
  fill_histos<<<blocks, p.threads>>>(d_an, d_events, d_results, p.process,
                                     n_events);
  sync_gpu_and_check("fill_histos");

  // copy the results back to the host
  cudaMemcpy(h_an, d_an, sizeof(analysis), cudaMemcpyDeviceToHost);

  // remove existing file
  std::remove(p.storage_file.c_str());

  // Calculate cross-section
  double xsec = h_an->wtot / h_an->ntot;      // Cross-section in pb
  double xsec_err = xsec / sqrt(h_an->ntot);  // Statistical error

  // Print out the total cross-section
  printf("Total cross-section: %.2e nb\n", xsec / 1000.);

  // Write cross-section in YODA format
  write_xsec(xsec, xsec_err, p.storage_file);

  // Scale and write histograms
  for (auto& hist : h_an->hists) {
    if (hist.name[0] != 'h') {
      // Special Case for L3: d59 is quoted per bin, not per unit n_ch, so the
      // bin width cancels the height = sumw / width that plotting applies
      if (std::string(hist.name) == "/L3_2004_I652683/d59-x01-y01") {
        hist.scale_w(hist.bins[0].width() / h_an->wtot);
      }

      // Other histos
      else {
        hist.scale_w(1. / h_an->wtot);
      }

      // Write to File
      write(hist, p.storage_file);
    }
  }

  // clean up
  delete h_an;
  cudaFree(d_an);
}