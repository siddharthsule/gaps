#include "observables.h"

// observable analysis

void analysis::analyze(event& ev) {
  // validate event
  ev.set_validity(ev.validate());

  if (!ev.get_validity()) {
    printf("invalid event\n");
    return;
  }

  // cluster
  cluster(ev);

  // calculate thrust
  calculate_thrust(ev);

  // calculate jet_m_br
  calculate_jet_m_br(ev);

  /**
   * why is the dalitz plot off?
   * ---------------------------
   *
   * while the dalitz analysis also benefits from the gpu parallelisation, the
   * writing of the data to file severely limits the performance, as instead of
   * the usual 100 bins, we have 100^2 = 1000 bins. this takes around 0.04s,
   * which is minute in the CPU case, but is in fact 40% of the total analysis
   * time! so for our tests, we keep this off, to keep our comparisons fair,
   * and relvant to the actual gpu effect.
   *
   * if you want to turn it on, uncomment the lines in this file, and it's
   * equivalent in the 'observables.cpp' file.
   */
  // calculate dalitz
  // calculate_dalitz(ev);

  // fill histograms
  hists[0].fill(ev.get_y23(), ev.get_dxs());
  hists[1].fill(ev.get_y34(), ev.get_dxs());
  hists[2].fill(ev.get_y45(), ev.get_dxs());
  hists[3].fill(ev.get_y56(), ev.get_dxs());
  hists[4].fill(ev.get_thr(), ev.get_dxs());
  hists[5].fill(ev.get_thr(), ev.get_dxs());
  hists[6].fill(ev.get_hjm(), ev.get_dxs());
  hists[7].fill(ev.get_ljm(), ev.get_dxs());
  hists[8].fill(ev.get_wjb(), ev.get_dxs());
  hists[9].fill(ev.get_njb(), ev.get_dxs());

  // dalitz plot is off
  // dalitz.fill(ev.get_dalitz(0), ev.get_dalitz(1), ev.get_dxs());

  // weighted total
  wtot += ev.get_dxs();
  ntot += 1.;
}

void analysis::finalize(const std::string& filename) {
  for (auto& hist : hists) {
    hist.scale_w(1. / ntot);
    hist.write(filename);
  }

  // dalitz plot is off
  // dalitz.scale_w(1. / ntot);
  // dalitz.write(filename);
}