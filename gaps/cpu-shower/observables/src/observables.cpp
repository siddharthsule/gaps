#include "observables.h"

void analysis::validate(event& ev) {
  /**
   * @brief Validate the event
   *
   * This function checks if the event is valid and sets the validity flag
   * accordingly.
   *
   * @param ev The event object
   */

  // validate event
  ev.set_validity(ev.validate());

  if (!ev.get_validity()) {
    printf("invalid event\n");
    return;
  }
}

void analysis::analyze(const event& ev) {
  /**
   * @brief Analyze the event and store the results
   *
   * This function first validates the event, then calculates the observables
   * for the corresponding process. The results are then stored in the
   * histograms and written to a file.
   *
   * @param ev The event object
   */

  // Make a vector to store the results of the analysis
  // For now, we set a size 20 per event and use that to store the results
  double results[21];
  for (int i = 0; i < 21; ++i) {
    results[i] = -50.0;
  }

  // jet rates
  cluster_durham(ev, results);

  // event shapes
  calculate_ev_shapes(ev, results);

  // fill histograms
  hists[0].fill(results[0], ev.get_dxs());
  hists[1].fill(results[1], ev.get_dxs());
  hists[2].fill(results[2], ev.get_dxs());
  hists[3].fill(results[3], ev.get_dxs());
  hists[4].fill(results[4], ev.get_dxs());
  hists[5].fill(results[5], ev.get_dxs());
  hists[6].fill(results[6], ev.get_dxs());
  hists[7].fill(results[7], ev.get_dxs());
  hists[8].fill(results[8], ev.get_dxs());
  hists[9].fill(results[9], ev.get_dxs());

  // weighted total
  wtot += ev.get_dxs();
  ntot += 1.;
}

void analysis::finalize(const std::string& filename) {
  /**
   * @brief Write the histograms to a file
   *
   * @param filename The name of the file to write the histograms to
   */

  for (auto& hist : hists) {
    if (hist.name != "hst") {
      hist.scale_w(1. / wtot);
      hist.write(filename);
    }
  }
}