#include "jetrates.h"

// Jet Rates

// Yij function Used for the Durham analysis
double Yij(const Vec4& p, const Vec4& q, double ecm2) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = std::min(p[0], q[0]);
  double max_pq = std::max(pq / std::sqrt(p.P2() * q.P2()), -1.);
  return 2. * std::pow(min_pq, 2) * (1. - std::min(max_pq, 1.)) / ecm2;
}

// Durham Clustering Algorithm
void Cluster(Event& ev) {
  if (!ev.GetValidity()) {
    return;
  }
  // For C++, one can use the std::vector class to store the parton 4-momenta.
  // Howver, this would be inefficient for CUDA, so we use a pre-set array size
  // To make the comparison fair, we make this class use the same array size
  // as well.

  // Get the center of mass energy squared
  double ecm2 = (ev.GetParton(0).GetMom() + ev.GetParton(1).GetMom()).M2();

  // Extract the 4-momenta of the partons
  Vec4 p[maxPartons];
  for (int i = 2; i < ev.GetSize(); ++i) {
    p[i - 2] = ev.GetParton(i).GetMom();
  }

  // kt2 will store the kt2 values for each clustering step
  // If not changed, set to -1 so we can ignore when histogramming
  double kt2[maxPartons] = {-1.};
  int counter = 0;

  // Number of partons (which will change when clustered), lower case to avoid N
  int n = ev.GetPartonSize();

  // imap will store the indices of the partons
  int imap[maxPartons];
  for (int i = 0; i < ev.GetPartonSize(); ++i) {
    imap[i] = i;
  }

  // kt2ij will store the kt2 values for each pair of partons
  double kt2ij[maxPartons][maxPartons] = {0.};
  double dmin = 1.;
  int ii = 0, jj = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double dij = kt2ij[i][j] = Yij(p[i], p[j], ecm2);
      if (dij < dmin) {
        dmin = dij;
        ii = i;
        jj = j;
      }
    }
  }

  // Cluster the partons
  while (n > 2) {
    --n;
    kt2[counter] = dmin;
    counter++;
    int jjx = imap[jj];
    p[jjx] = p[jjx] + p[imap[ii]];
    for (int i = ii; i < n; ++i) {
      imap[i] = imap[i + 1];
    }
    for (int j = 0; j < jj; ++j) {
      kt2ij[jjx][imap[j]] = Yij(p[jjx], p[imap[j]], ecm2);
    }
    for (int i = jj + 1; i < n; ++i) {
      kt2ij[imap[i]][jjx] = Yij(p[jjx], p[imap[i]], ecm2);
    }
    dmin = 1.;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < i; ++j) {
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
  ev.SetY23(counter > 0 ? log10(kt2[counter - 1 - 0]) : -50.);
  ev.SetY34(counter > 1 ? log10(kt2[counter - 1 - 1]) : -50.);
  ev.SetY45(counter > 2 ? log10(kt2[counter - 1 - 2]) : -50.);
  ev.SetY56(counter > 3 ? log10(kt2[counter - 1 - 3]) : -50.);
}
