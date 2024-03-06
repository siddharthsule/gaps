#include "observables.h"

// -----------------------------------------------------------------------------
// Jet Rates

double Yij(const Vec4& p, const Vec4& q, double ecm2) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = std::min(p[0], q[0]);
  double max_pq = std::max(pq / std::sqrt(p.P2() * q.P2()), -1.0);
  return 2.0 * std::pow(min_pq, 2) * (1.0 - std::min(max_pq, 1.0)) / ecm2;
}

void Cluster(Event& ev) {
  double ecm2 = (ev.GetParton(0).GetMom() + ev.GetParton(1).GetMom()).M2();
  std::vector<Vec4> p;
  for (size_t i = 2; i < ev.GetSize(); ++i) {
    p.push_back(ev.GetParton(i).GetMom());
  }

  std::vector<double> kt2;
  size_t N = p.size();
  std::vector<size_t> imap(N);
  for (size_t i = 0; i < N; ++i) {
    imap[i] = i;
  }

  std::vector<std::vector<double>> kt2ij(N, std::vector<double>(N, 0.0));
  double dmin = 1.0;
  size_t ii = 0, jj = 0;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < i; ++j) {
      double dij = kt2ij[i][j] = Yij(p[i], p[j], ecm2);
      if (dij < dmin) {
        dmin = dij;
        ii = i;
        jj = j;
      }
    }
  }

  while (N > 2) {
    --N;
    kt2.push_back(dmin);
    size_t jjx = imap[jj];
    p[jjx] = p[jjx] + p[imap[ii]];
    for (size_t i = ii; i < N; ++i) {
      imap[i] = imap[i + 1];
    }
    for (size_t j = 0; j < jj; ++j) {
      kt2ij[jjx][imap[j]] = Yij(p[jjx], p[imap[j]], ecm2);
    }
    for (size_t i = jj + 1; i < N; ++i) {
      kt2ij[imap[i]][jjx] = Yij(p[jjx], p[imap[i]], ecm2);
    }
    dmin = 1.0;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < i; ++j) {
        double dij = kt2ij[imap[i]][imap[j]];
        if (dij < dmin) {
          dmin = dij;
          ii = i;
          jj = j;
        }
      }
    }
  }

  // Store the kt2 values in the output arrays
  if (!kt2.empty()) {
    ev.SetY23(kt2.size() > 0 ? std::log10(kt2[kt2.size() - 1 - 0]) : -50.0);
    ev.SetY34(kt2.size() > 1 ? std::log10(kt2[kt2.size() - 1 - 1]) : -50.0);
    ev.SetY45(kt2.size() > 2 ? std::log10(kt2[kt2.size() - 1 - 2]) : -50.0);
    ev.SetY56(kt2.size() > 3 ? std::log10(kt2[kt2.size() - 1 - 3]) : -50.0);
  }
}

// ----------------------------------------------------------------------------- 
// Event Shapes

std::vector<Vec4> GetMomenta(const Event& ev) {
  std::vector<Vec4> moms;

  // Get Momenta for partons starting at 2
  for (size_t i = 2; i < ev.GetSize(); ++i) {
    Vec4 pmom = ev.GetParton(i).GetMom();
    moms.push_back(pmom);
  }

  // Sort by mag2 in descending order (of 3 Momenta)
  std::sort(moms.begin(), moms.end(),
            [](const Vec4& a, const Vec4& b) { return a.P() > b.P(); });

  return moms;
}

void CalculateThrust(Event& ev) {
  std::vector<Vec4> moms = GetMomenta(ev);

  if (moms.size() < 3) {
    return;
  }

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

  ev.SetThr(thrust);
  ev.SetTAxis(t_axis);
}

void CalculateJetMBr(Event& ev) {
  std::vector<Vec4> moms = GetMomenta(ev);

  if (moms.size() < 3) {
    return;
  }

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
    double mo_para = mom.Dot(ev.GetTAxis());
    double mo_perp = (mom - (ev.GetTAxis() * mo_para)).P();
    double enrg = mom.P();  // Equivalent to mom.E for massless partons

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
      p_against = p_against + (mom * 0.5);
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
    ev.SetHJM(mH);
    ev.SetWJB(bW);
  } else {
    ev.SetHJM(mH);
    ev.SetLJM(mL);
    ev.SetWJB(bW);
    ev.SetNJB(bN);
  }
}

// -----------------------------------------------------------------------------
// Observable Analysis

void Analysis::Analyze(Event& ev) {

  // Cluster
  Cluster(ev);

  // Calculate Thrust
  CalculateThrust(ev);

  // Calculate JetMBr
  CalculateJetMBr(ev);

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

  // Weighted Total
  wtot += ev.GetDxs();
  ntot += 1.0;
}

void Analysis::Finalize(const std::string& filename) {
  for (auto& hist : hists) {
    hist.ScaleW(1.0 / ntot);
    hist.Write(filename);
  }
}