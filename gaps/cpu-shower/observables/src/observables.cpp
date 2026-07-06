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
  double results[20];
  for (int i = 0; i < 20; ++i) {
    results[i] = -50.0;
  }

  // LEP: e+ e- -> q qbar
  if (process == 1) {
    // jet rates
    double jet_rates[4] = {-50., -50., -50., -50.};
    cluster_durham(ev, jet_rates);
    double log10y23 = jet_rates[0];
    double log10y34 = jet_rates[1];
    double log10y45 = jet_rates[2];
    double log10y56 = jet_rates[3];

    // event shapes
    double ev_shapes[5] = {-50., -50., -50., -50., -50.};
    calculate_ev_shapes(ev, ev_shapes);
    double thrust = ev_shapes[0];
    double hjm = ev_shapes[1];
    double ljm = ev_shapes[2];
    double wjb = ev_shapes[3];
    double njb = ev_shapes[4];

    // fill histograms
    hists[0].fill(log10y23, ev.get_dxs());
    hists[1].fill(log10y34, ev.get_dxs());
    hists[2].fill(log10y45, ev.get_dxs());
    hists[3].fill(log10y56, ev.get_dxs());
    hists[4].fill(1. - thrust, ev.get_dxs());
    hists[5].fill(1. - thrust, ev.get_dxs());
    hists[6].fill(hjm, ev.get_dxs());
    hists[7].fill(ljm, ev.get_dxs());
    hists[8].fill(wjb, ev.get_dxs());
    hists[9].fill(njb, ev.get_dxs());
    hists[10].fill(ev.get_size() - 2, ev.get_dxs());

    // ALEPH
    hists[11].fill(thrust, ev.get_dxs());
    hists[12].fill(sqr(hjm), ev.get_dxs());
    hists[13].fill(wjb, ev.get_dxs());
    hists[14].fill(sqr(hjm) - sqr(ljm), ev.get_dxs());
    hists[15].fill(wjb + njb, ev.get_dxs());
    hists[16].fill(-log(pow(10, log10y23)), ev.get_dxs());
    hists[17].fill(-log(pow(10, log10y34)), ev.get_dxs());
    hists[18].fill(-log(pow(10, log10y45)), ev.get_dxs());
    hists[19].fill(-log(pow(10, log10y56)), ev.get_dxs());

    // L3
    double cmul[1] = {-50.};
    calculate_chargedmult(ev, cmul);
    hists[20].fill(cmul[0], ev.get_dxs());
  }

  // LHC: p p -> e+ e-
  else if (process == 2) {
    // calculate Z boson observables
    calculate_mczinc(ev, results);

    // Gen kt Algorithm
    cluster_genkt(ev, results);

    // fill histograms
    hists[0].fill(results[0], ev.get_dxs());
    hists[1].fill(results[1], ev.get_dxs());
    hists[2].fill(results[1], ev.get_dxs());
    hists[3].fill(results[2], ev.get_dxs());
    hists[4].fill(results[3], ev.get_dxs());
    hists[5].fill(results[4], ev.get_dxs() / 2);
    hists[5].fill(results[5], ev.get_dxs() / 2);
    hists[6].fill(results[6], ev.get_dxs() / 2);
    hists[6].fill(results[7], ev.get_dxs() / 2);
    hists[7].fill(results[8], ev.get_dxs());
    hists[8].fill(results[9], ev.get_dxs());
    hists[9].fill(results[10], ev.get_dxs());
    hists[10].fill(results[11], ev.get_dxs());
    hists[11].fill(results[12], ev.get_dxs());
    hists[12].fill(results[13], ev.get_dxs());
    hists[13].fill(results[14], ev.get_dxs());
    hists[14].fill(results[15], ev.get_dxs());
    hists[15].fill(results[16], ev.get_dxs());
    hists[16].fill(results[17], ev.get_dxs());
    hists[17].fill(results[18], ev.get_dxs());
    hists[18].fill(results[19], ev.get_dxs());

    // Forward-Backward Asymmetry testing
    // hists[3].fill(-results[3], -ev.get_dxs());
  }

  // weighted total
  wtot += ev.get_dxs();
  ntot += 1.;
}

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

void analysis::finalize(const std::string& filename) {
  /**
   * @brief Write the histograms to a file
   *
   * @param filename The name of the file to write the histograms to
   */

  // Calculate cross-section
  double xsec = wtot / ntot;            // Cross-section in pb
  double xsec_err = xsec / sqrt(ntot);  // Statistical error

  // Print out the total cross-section
  printf("Total cross-section: %.2e nb\n", xsec / 1000.);

  // Write cross-section in YODA format
  write_xsec(xsec, xsec_err, filename);

  // Scale and write histograms
  for (auto& hist : hists) {
    if (hist.name != "hst") {
      hist.scale_w(1. / wtot);
      hist.write(filename);
    }
  }
}