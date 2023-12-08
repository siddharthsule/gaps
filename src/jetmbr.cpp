#include "jetmbr.h"
#include <tuple>

JetMBr::JetMBr() {}

std::vector<double> JetMBr::CalculateJetMBr(const Event& ev) {
  std::vector<double> jetmbr = {-5.0, -5.0, -5.0, -5.0};
  std::vector<std::vector<double>> moms = GetMomenta(ev);
  std::tuple<double, std::vector<double>> thrust = thr.CalculateThrust(ev);
  std::vector<double> t_axis = std::get<1>(thrust);

  // Momentum Sum (Should be 91.2 GeV)
  double momsum = 0.0;
  for (const auto& mom : moms) {
    momsum += Magnitude(mom);
  }

  Vec4 p_with, p_against;
  int n_with = 0, n_against = 0;
  double e_vis = 0.0, broad_with = 0.0, broad_against = 0.0, broad_denominator = 0.0;

  for (const auto& mom : moms) {
    double mo_para = DotProduct(mom, t_axis);
    double mo_perp = Magnitude(Subtract(mom, Multiply(t_axis, mo_para)));
    double enrg = Magnitude(mom);

    e_vis += enrg;
    broad_denominator += 2.0 * enrg;

    Vec4 temp = Vec4(Magnitude(mom), mom[0], mom[1], mom[2]);

    if (mo_para > 0.0) {
      p_with = p_with + temp;
      broad_with += mo_perp;
      n_with++;
    } else if (mo_para < 0.0) {
      p_against = p_against + temp;
      broad_against += mo_perp;
      n_against++;
    } else {
      p_with = p_with + (Vec4(Magnitude(mom), mom[0], mom[1], mom[2]) * 0.5);
      p_against = p_against + (Vec4(Magnitude(mom), mom[0], mom[1], mom[2]) * 0.5);
      broad_with += 0.5 * mo_perp;
      broad_against += 0.5 * mo_perp;
      n_with++;
      n_against++;
    }
  }

  double e2_vis = std::pow(e_vis, 2.0);

  double mass2_with = std::abs(p_with.M2() / e2_vis);
  double mass2_against = std::abs(p_against.M2() / e2_vis);

  double mass_with = std::sqrt(mass2_with);
  double mass_against = std::sqrt(mass2_against);

  broad_with /= broad_denominator;
  broad_against /= broad_denominator;

  double mH = std::max(mass_with, mass_against);
  double mL = std::min(mass_with, mass_against);

  double bW = std::max(broad_with, broad_against);
  double bN = std::min(broad_with, broad_against);

  if (n_with == 1 || n_against == 1) {
    jetmbr[0] = mH;
    jetmbr[1] = -5.0;
    jetmbr[2] = bW;
    jetmbr[3] = -5.0;
  } else {
    jetmbr[0] = mH;
    jetmbr[1] = mL;
    jetmbr[2] = bW;
    jetmbr[3] = bN;
  }

  return jetmbr;
}

MBrAnalysis::MBrAnalysis() : jetmbr(), wtot(0.0) {
  jetmbr_hists.push_back(Histo1D(100, 0.0, 1.0, "/GPU_EvGEN/hjm"));
  jetmbr_hists.push_back(Histo1D(100, 0.0, 0.5, "/GPU_EvGEN/ljm"));
  jetmbr_hists.push_back(Histo1D(100, 0.0, 0.5, "/GPU_EvGEN/wjb"));
  jetmbr_hists.push_back(Histo1D(100, 0.0, 0.2, "/GPU_EvGEN/njb"));
}

void MBrAnalysis::Analyze(const Event& ev) {
  std::vector<double> mbrdata = jetmbr.CalculateJetMBr(ev);

  jetmbr_hists[0].Fill(mbrdata[0], ev.GetDxs());
  jetmbr_hists[1].Fill(mbrdata[1], ev.GetDxs());
  jetmbr_hists[2].Fill(mbrdata[2], ev.GetDxs());
  jetmbr_hists[3].Fill(mbrdata[3], ev.GetDxs());

  wtot += ev.GetDxs();
}

void MBrAnalysis::Finalize(const std::string& filename) {
  for (auto& hist : jetmbr_hists) {
    hist.ScaleW(1.0 / wtot);
    hist.Write(filename);
  }
}