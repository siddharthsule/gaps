#include "observables.h"

// Observable Analysis

void Analysis::Analyze(Event& ev) {
  // Validate Event
  ev.SetValidity(ev.Validate());

  if (!ev.GetValidity()) {
    printf("Invalid Event\n");
    return;
  }

  // Cluster
  Cluster(ev);

  // Calculate Thrust
  CalculateThrust(ev);

  // Calculate JetMBr
  CalculateJetMBr(ev);

  // Calculate Dalitz
  CalculateDalitz(ev);

  // Fill Histograms
  hists[0].Fill(ev.GetY23(), ev.GetDxs());
  hists[1].Fill(ev.GetY34(), ev.GetDxs());
  hists[2].Fill(ev.GetY45(), ev.GetDxs());
  hists[3].Fill(ev.GetY56(), ev.GetDxs());
  hists[4].Fill(ev.GetThr(), ev.GetDxs());
  hists[5].Fill(ev.GetThr(), ev.GetDxs());
  hists[6].Fill(ev.GetHJM(), ev.GetDxs());
  hists[7].Fill(ev.GetLJM(), ev.GetDxs());
  hists[8].Fill(ev.GetWJB(), ev.GetDxs());
  hists[9].Fill(ev.GetNJB(), ev.GetDxs());

  dalitz.Fill(ev.GetDalitz(0), ev.GetDalitz(1), ev.GetDxs());

  // Weighted Total
  wtot += ev.GetDxs();
  ntot += 1.;
}

void Analysis::Finalize(const std::string& filename) {
  for (auto& hist : hists) {
    hist.ScaleW(1. / ntot);
    hist.Write(filename);
  }
  dalitz.ScaleW(1. / ntot);
  dalitz.Write(filename);
}