#include <fstream>

#include "observables.cuh"

// -----------------------------------------------------------------------------
// validate events before binning

__global__ void validate_events(event* events, int* invalid, int n) {
  int idx = block_idx.x * block_dim.x + thread_idx.x;

  if (idx >= n) {
    return;
  }

  event& ev = events[idx];
  ev.set_validity(ev.validate());

  if (!ev.get_validity()) {
    // printf("invalid event\n");
    atomic_add(invalid, 1);
  }
}

// -----------------------------------------------------------------------------
// analysis

// fill the_histograms (atomically!)
__global__ void fill_histos(analysis* an, event* events, int n) {
  int idx = block_idx.x * block_dim.x + thread_idx.x;

  if (idx >= n) {
    return;
  }

  event& ev = events[idx];

  an->hists[0].fill(ev.get_y23(), ev.get_dxs());
  an->hists[1].fill(ev.get_y34(), ev.get_dxs());
  an->hists[2].fill(ev.get_y45(), ev.get_dxs());
  an->hists[3].fill(ev.get_y56(), ev.get_dxs());
  an->hists[4].fill(ev.get_thr(), ev.get_dxs());
  an->hists[5].fill(ev.get_thr(), ev.get_dxs());
  an->hists[6].fill(ev.get_hjm(), ev.get_dxs());
  an->hists[7].fill(ev.get_ljm(), ev.get_dxs());
  an->hists[8].fill(ev.get_wjb(), ev.get_dxs());
  an->hists[9].fill(ev.get_njb(), ev.get_dxs());

  // dalitz plot is off
  // an->dalitz.fill(ev.get_dalitz(0), ev.get_dalitz(1), ev.get_dxs());

  atomic_add(&an->wtot, ev.get_dxs());
  atomic_add(&an->ntot, 1.);
}

// run the above kernels
void do_analysis(thrust::device_vector<event>& d_events, std::string filename) {
  /**
   * only place for a host object
   * ----------------------------
   *
   * while most of the work is done on the device, one cannot directly write
   * to a file from the device. therefore, we will create a host object to
   * store the histograms, and then copy the results back to the host for
   * writing to file.
   */

  // device analysis object
  analysis *h_an, *d_an;

  // allocate memory for the device analysis object
  h_an = new analysis();
  cuda_malloc(&d_an, sizeof(analysis));
  cuda_memcpy(d_an, h_an, sizeof(analysis), cuda_memcpy_host_to_device);

  // get event data
  int n = d_events.size();
  event* d_events_ptr = thrust::raw_pointer_cast(d_events.data());

  // validate the events
  int* d_invalid;
  cuda_malloc(&d_invalid, sizeof(int));
  cuda_memset(d_invalid, 0, sizeof(int));

  validate_events<<<(n + 255) / 256, 256>>>(d_events_ptr, d_invalid, n);
  sync_gpu_and_check("validate_events");

  int h_invalid;
  cuda_memcpy(&h_invalid, d_invalid, sizeof(int), cuda_memcpy_device_to_host);
  cuda_free(d_invalid);

  if (h_invalid > 0) {
    std::cout << "" << std::endl;
    std::cout << "error: invalid events found" << std::endl;
    std::cout << "number of invalid events: " << h_invalid << "\n";
  }

  // calculare the observables
  do_cluster<<<(n + 255) / 256, 256>>>(d_events_ptr, n);
  sync_gpu_and_check("do_cluster");

  calculate_thr<<<(n + 255) / 256, 256>>>(d_events_ptr, n);
  sync_gpu_and_check("calculate_thr");

  calculate_jet_m_br<<<(n + 255) / 256, 256>>>(d_events_ptr, n);
  sync_gpu_and_check("calculate_jet_m_br");

  /**
   * why is the dalitz plot off?
   * ---------------------------
   *
   * while the dalitz analysis also benefits from the gpu parallelisation, the
   * writing of the data to file severely limits the performance, as instead of
   * the usual 100 bins, we have 100^2 = 1000 bins. this takes around 0.04s,
   * which is minute in the c++ case, but is in fact 40% of the total analysis
   * time! so for our tests, we keep this off, to keep our comparisons fair,
   * and relvant to the actual gpu effect.
   *
   * if you want to turn it on, uncomment the lines in this file, and it's
   * equivalent in the 'observables.cpp' file.
   */
  // calculate_dalitz<<<(n + 255) / 256, 256>>>(d_events_ptr, n);
  // sync_gpu_and_check("calculate_dalitz");

  // do the analysis
  fill_histos<<<(n + 255) / 256, 256>>>(d_an, d_events_ptr, n);
  sync_gpu_and_check("fill_histos");

  // copy the results back to the host
  cuda_memcpy(h_an, d_an, sizeof(analysis), cuda_memcpy_device_to_host);

  // normalize the histograms
  for (auto& hist : h_an->hists) {
    hist.scale_w(1. / h_an->ntot);
  }

  // dalitz plot is off
  // h_an->dalitz.scale_w(1. / h_an->ntot);

  // remove existing file
  std::remove(filename.c_str());

  // write the histograms to file
  write(h_an->hists[0], "/gaps/log10y23\n", filename);
  write(h_an->hists[1], "/gaps/log10y34\n", filename);
  write(h_an->hists[2], "/gaps/log10y45\n", filename);
  write(h_an->hists[3], "/gaps/log10y56\n", filename);
  write(h_an->hists[4], "/gaps/tvalue\n", filename);
  write(h_an->hists[5], "/gaps/tzoomd\n", filename);
  write(h_an->hists[6], "/gaps/hjm\n", filename);
  write(h_an->hists[7], "/gaps/ljm\n", filename);
  write(h_an->hists[8], "/gaps/wjb\n", filename);
  write(h_an->hists[9], "/gaps/njb\n", filename);

  // dalitz plot is off
  // write(h_an->dalitz, "/gaps/dalitz\n", filename);

  // clean up
  delete h_an;
  cuda_free(d_an);
}