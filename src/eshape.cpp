#include "eshape.h"

#include <algorithm>
#include <cmath>
#include <vector>

EShape::EShape() {}

std::vector<Vec4> EShape::GetMomenta(const Event& ev) {
  std::vector<Vec4> moms;

  // Get Momenta for particles starting at 2
  for (size_t i = 2; i < ev.GetSize(); ++i) {
    Vec4 pmom = ev.GetParton(i).GetMom();
    moms.push_back(pmom);
  }

  // Sort by mag2 in descending order (of 3 Momenta)
  std::sort(moms.begin(), moms.end(),
            [](const Vec4& a, const Vec4& b) { return a.P() > b.P(); });

  return moms;
}

std::tuple<double, Vec4> EShape::CalculateThrust(const Event& ev) {
  std::vector<Vec4> moms = GetMomenta(ev);

  // Sum of all momenta
  double momsum = 0.0;
  for (auto& mom : moms) {
    momsum += mom.P();
  }

  double thrust = 0.0;
  Vec4 t_axis = Vec4();

  for (size_t k = 1; k < moms.size(); ++k) {
    for (size_t j = 0; j < k; ++j) {
      Vec4 tmp_axis = moms[j].Cross(moms[k]);
      Vec4 p_thrust = Vec4();
      std::vector<Vec4> p_combin;

      for (size_t i = 0; i < moms.size(); ++i) {
        if (i != j && i != k) {
          if (moms[i].Dot(tmp_axis) >= 0) {
            p_thrust = p_thrust + moms[i];
          } else {
            p_thrust = p_thrust - moms[i];
          }
        }
      }

      p_combin.push_back((p_thrust + moms[j] + moms[k]));
      p_combin.push_back((p_thrust + moms[j] - moms[k]));
      p_combin.push_back((p_thrust - moms[j] + moms[k]));
      p_combin.push_back((p_thrust - moms[j] - moms[k]));

      for (auto p : p_combin) {
        double temp = p.P();
        if (temp > thrust) {
          thrust = temp;
          t_axis = p;  // will unit-ify below
        }
      }
    }
  }

  thrust /= momsum;
  thrust = 1.0 - thrust;

  t_axis = t_axis / (t_axis).P();
  if (t_axis[2] < 0) {
    t_axis = t_axis * -1.0;
  }

  if (thrust < 1e-12) {
    thrust = -5.0;
  }

  return std::make_tuple(thrust, t_axis);
}

std::vector<double> EShape::CalculateJetMBr(const Event& ev) {
  
  
  std::vector<double> jetmbr = {-5.0, -5.0, -5.0, -5.0};
  std::vector<Vec4> moms = GetMomenta(ev);
  std::tuple<double, Vec4> thrust = CalculateThrust(ev);
  Vec4 t_axis = std::get<1>(thrust);

  // Momentum Sum (Should be 91.2 GeV)
  double momsum = 0.0;
  for (const auto& mom : moms) {
    momsum += mom.P();
  }

  Vec4 p_with, p_against;
  int n_with = 0, n_against = 0;
  double e_vis = 0.0, broad_with = 0.0, broad_against = 0.0,
         broad_denominator = 0.0;

  for (const auto& mom : moms) {
    double mo_para = mom.Dot(t_axis);
    double mo_perp = (mom - (t_axis * mo_para)).P();
    double enrg = mom.P(); // Equivalent to mom.E for massless particles

    e_vis += enrg;
    broad_denominator += 2.0 * enrg;

    if (mo_para > 0.0) {
      p_with = p_with + mom;
      broad_with += mo_perp;
      n_with++;
    } else if (mo_para < 0.0) {
      p_against = p_against + mom;
      broad_against += mo_perp;
      n_against++;
    } else {
      p_with = p_with + (mom * 0.5);
      p_against =
          p_against + (mom * 0.5);
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

EAnalysis::EAnalysis() : eshape(), wtot(0.0) {
  eshape_hists.push_back(Histo1D(100, 0.0, 0.5, "/GPU_EvGEN/tvalue"));
  eshape_hists.push_back(Histo1D(100, 0.0, 0.5, "/GPU_EvGEN/tzoomd"));
  eshape_hists.push_back(Histo1D(100, 0.0, 1.0, "/GPU_EvGEN/hjm"));
  eshape_hists.push_back(Histo1D(100, 0.0, 0.5, "/GPU_EvGEN/ljm"));
  eshape_hists.push_back(Histo1D(100, 0.0, 0.5, "/GPU_EvGEN/wjb"));
  eshape_hists.push_back(Histo1D(100, 0.0, 0.2, "/GPU_EvGEN/njb"));
}

void EAnalysis::Analyze(const Event& ev) {
  std::tuple<double, Vec4> thrust = eshape.CalculateThrust(ev);
  std::vector<double> mbrdata = eshape.CalculateJetMBr(ev);

  eshape_hists[0].Fill(std::get<0>(thrust), ev.GetDxs());
  eshape_hists[1].Fill(std::get<0>(thrust), ev.GetDxs());
  eshape_hists[2].Fill(mbrdata[0], ev.GetDxs());
  eshape_hists[3].Fill(mbrdata[1], ev.GetDxs());
  eshape_hists[4].Fill(mbrdata[2], ev.GetDxs());
  eshape_hists[5].Fill(mbrdata[3], ev.GetDxs());

  wtot += ev.GetDxs();
}

void EAnalysis::Finalize(const std::string& filename) {
  for (auto& hist : eshape_hists) {
    hist.ScaleW(1.0 / wtot);
    hist.Write(filename);
  }
}