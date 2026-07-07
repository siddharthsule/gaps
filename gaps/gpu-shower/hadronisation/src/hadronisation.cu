#include "hadronisation.cuh"

void run_hadronisation(thrust::device_vector<event>& dv_events, const params& p,
                       int blocks) {
  /**
   * @brief Run the hadronisation on the events
   *
   * Runs the five hadronisation steps in sequence, one kernel per step.
   * Each kernel processes all events in parallel (one thread per event).
   *
   * @param dv_events The device vector of showered events
   * @param p         Run parameters (threads, clmax, clpow, psplit, ...)
   * @param blocks    Number of CUDA blocks
   */

  // Get raw pointer and event count from the thrust vector
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // Set up the device hadronisation object
  hadronisation *h_had, *d_had;
  h_had = new hadronisation(p.clmax, p.clpow, p.psplit, p.pwt);
  cudaMalloc(&d_had, sizeof(hadronisation));
  cudaMemcpy(d_had, h_had, sizeof(hadronisation), cudaMemcpyHostToDevice);

  // Allocate the cluster lists — one per event, only live during hadronisation
  thrust::device_vector<cluster_list> dv_cls(n_events);
  cluster_list* d_cls = thrust::raw_pointer_cast(dv_cls.data());

  // Run the hadronisation steps in sequence
  debug_msg("running @constituent_reshuffling");
  constituent_reshuffling<<<blocks, p.threads>>>(d_events, d_had, n_events);
  sync_gpu_and_check("constituent_reshuffling");

  debug_msg("running @force_gluons_to_split");
  force_gluons_to_split<<<blocks, p.threads>>>(d_events, d_had, n_events);
  sync_gpu_and_check("force_gluons_to_split");

  debug_msg("running @form_clusters");
  form_clusters<<<blocks, p.threads>>>(d_events, d_cls, d_had, n_events);
  sync_gpu_and_check("form_clusters");

  debug_msg("running @fission_clusters");
  fission_clusters<<<blocks, p.threads>>>(d_events, d_cls, d_had, n_events);
  sync_gpu_and_check("fission_clusters");

  debug_msg("running @decay_clusters");
  decay_clusters<<<blocks, p.threads>>>(d_events, d_cls, d_had, n_events);
  sync_gpu_and_check("decay_clusters");

  // Clean up (dv_cls is freed automatically when it goes out of scope)
  delete h_had;
  cudaFree(d_had);
}
