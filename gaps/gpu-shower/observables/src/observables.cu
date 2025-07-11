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
                            int n) {
  /**
   * @brief Fill the histograms with the results of the analysis
   *
   * @param an The analysis object
   * @param ev The event object
   * @param results The array to store the results
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------

  const event& ev = events[idx];

  // I know forloop can be done here, but just to note what is happening
  an->hists[0].fill(results[21 * idx + 0], ev.get_dxs());  // log10y23
  an->hists[1].fill(results[21 * idx + 1], ev.get_dxs());  // log10y34
  an->hists[2].fill(results[21 * idx + 2], ev.get_dxs());  // log10y45
  an->hists[3].fill(results[21 * idx + 3], ev.get_dxs());  // log10y56
  an->hists[4].fill(results[21 * idx + 4], ev.get_dxs());  // tvalue
  an->hists[5].fill(results[21 * idx + 5], ev.get_dxs());  // tzoomd
  an->hists[6].fill(results[21 * idx + 6], ev.get_dxs());  // hjm
  an->hists[7].fill(results[21 * idx + 7], ev.get_dxs());  // ljm
  an->hists[8].fill(results[21 * idx + 8], ev.get_dxs());  // wjb
  an->hists[9].fill(results[21 * idx + 9], ev.get_dxs());  // njb

  atomicAdd(&an->wtot, ev.get_dxs());
  atomicAdd(&an->ntot, 1.);
}

// run the above kernels
void do_analysis(thrust::device_vector<event>& dv_events, std::string filename,
                 int blocks, int threads) {
  /**
   * @brief Run the analysis
   *
   * @param dv_events device vector of event records
   * @param filename output file name
   * @param blocks number of blocks to use for the kernel
   * @param threads number of threads per block
   */

  // device analysis object
  analysis *h_an, *d_an;

  // allocate memory for the device analysis object
  h_an = new analysis();
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

  // do the analysis

  // Make a vector to store the results of the analysis
  // For now, we set a size 10 per event and use that to store the results
  thrust::device_vector<double> dv_results(21 * n_events, -50.);
  double* d_results = thrust::raw_pointer_cast(dv_results.data());

  cluster_durham<<<blocks, threads>>>(d_events, d_results, n_events);
  sync_gpu_and_check("cluster_durham");

  calculate_ev_shapes<<<blocks, threads>>>(d_events, d_results, n_events);
  sync_gpu_and_check("calculate_ev_shapes");

  // do the analysis
  fill_histos<<<blocks, threads>>>(d_an, d_events, d_results, n_events);
  sync_gpu_and_check("fill_histos");

  // copy the results back to the host
  cudaMemcpy(h_an, d_an, sizeof(analysis), cudaMemcpyDeviceToHost);

  // normalize the histograms
  for (auto& hist : h_an->hists) {
    hist.scale_w(1. / h_an->wtot);
  }

  // remove existing file
  std::remove(filename.c_str());

  for (auto& hist : h_an->hists) {
    if (hist.name[0] != 'h') {
      write(hist, filename);
    }
  }

  // clean up
  delete h_an;
  cudaFree(d_an);
}