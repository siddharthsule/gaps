#include "durham.h"

#include <cmath>

Durham::Durham() {}

double Durham::Yij(const Vec4& p, const Vec4& q) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = std::min(p[0], q[0]);
  double max_pq = std::max(pq / std::sqrt(p.P2() * q.P2()), -1.0);
  return 2.0 * std::pow(min_pq, 2) * (1.0 - std::min(max_pq, 1.0)) / ecm2;
}

std::vector<double> Durham::Cluster(const Event& ev) {
  ecm2 = (ev.GetParton(0).GetMom() + ev.GetParton(1).GetMom()).M2();
  std::vector<Vec4> p;
  for (size_t i = 2; i < ev.GetSize(); ++i) {
    p.push_back(ev.GetParton(i).GetMom());
  }

  std::vector<double> kt2;
  size_t n = p.size();
  std::vector<size_t> imap(n);
  for (size_t i = 0; i < n; ++i) {
    imap[i] = i;
  }

  std::vector<std::vector<double>> kt2ij(n, std::vector<double>(n, 0.0));
  double dmin = 1.0;
  size_t ii = 0, jj = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      double dij = kt2ij[i][j] = Yij(p[i], p[j]);
      if (dij < dmin) {
        dmin = dij;
        ii = i;
        jj = j;
      }
    }
  }

  while (n > 2) {
    --n;
    kt2.push_back(dmin);
    size_t jjx = imap[jj];
    p[jjx] = p[jjx] + p[imap[ii]];
    for (size_t i = ii; i < n; ++i) {
      imap[i] = imap[i + 1];
    }
    for (size_t j = 0; j < jj; ++j) {
      kt2ij[jjx][imap[j]] = Yij(p[jjx], p[imap[j]]);
    }
    for (size_t i = jj + 1; i < n; ++i) {
      kt2ij[imap[i]][jjx] = Yij(p[jjx], p[imap[i]]);
    }
    dmin = 1.0;
    for (size_t i = 0; i < n; ++i) {
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
  return kt2;
}

DAnalysis::DAnalysis() : dur(), wtot(0.0) {
  dur_hists.emplace_back(100, -5.0, 0.0, "/GPU_EvGEN/log10y23\n");
  dur_hists.emplace_back(100, -5.0, 0.0, "/GPU_EvGEN/log10y34\n");
  dur_hists.emplace_back(100, -5.0, 0.0, "/GPU_EvGEN/log10y45\n");
  dur_hists.emplace_back(100, -5.0, 0.0, "/GPU_EvGEN/log10y56\n");
}

void DAnalysis::Analyze(const Event& ev) {
  std::vector<double> kt2 = dur.Cluster(ev);

  for (size_t j = 0; j < dur_hists.size(); ++j) {
    dur_hists[j].Fill(
        (kt2.size() > j) ? std::log10(kt2[kt2.size() - 1 - j]) : -50.0,
        ev.GetDxs());
  }

  wtot += ev.GetDxs();
}

void DAnalysis::Finalize(const std::string& filename) {
  for (auto& dur_hist : dur_hists) {
    dur_hist.ScaleW(1.0 / wtot);
    dur_hist.Write(filename);
  }
}