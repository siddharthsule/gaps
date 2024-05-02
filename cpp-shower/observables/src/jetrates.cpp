#include "jetrates.h"

// jet rates

// yij function used for the durham analysis
double yij(const vec4& p, const vec4& q, double ecm2) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = std::min(p[0], q[0]);
  double max_pq = std::max(pq / std::sqrt(p.p2() * q.p2()), -1.);
  return 2. * std::pow(min_pq, 2) * (1. - std::min(max_pq, 1.)) / ecm2;
}

// durham clustering algorithm
void cluster(event& ev) {
  if (!ev.get_validity()) {
    return;
  }
  // for c++, one can use the std::vector class to store the parton 4-momenta.
  // howver, this would be inefficient for cuda, so we use a pre-set array size
  // to make the comparison fair, we make this class use the same array size
  // as well.

  // get the center of mass energy squared
  double ecm2 = (ev.get_parton(0).get_mom() + ev.get_parton(1).get_mom()).m2();

  // extract the 4-momenta of the partons
  vec4 p[max_partons];
  for (int i = 2; i < ev.get_size(); ++i) {
    p[i - 2] = ev.get_parton(i).get_mom();
  }

  // kt2 will store the kt2 values for each clustering step
  // if not changed, set to -1 so we can ignore when histogramming
  double kt2[max_partons] = {-1.};
  int counter = 0;

  // number of partons (which will change when clustered), lower case to avoid n
  int n = ev.get_parton_size();

  // imap will store the indices of the partons
  int imap[max_partons];
  for (int i = 0; i < ev.get_parton_size(); ++i) {
    imap[i] = i;
  }

  // kt2ij will store the kt2 values for each pair of partons
  double kt2ij[max_partons][max_partons] = {0.};
  double dmin = 1.;
  int ii = 0, jj = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      double dij = kt2ij[i][j] = yij(p[i], p[j], ecm2);
      if (dij < dmin) {
        dmin = dij;
        ii = i;
        jj = j;
      }
    }
  }

  // cluster the partons
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
      kt2ij[jjx][imap[j]] = yij(p[jjx], p[imap[j]], ecm2);
    }
    for (int i = jj + 1; i < n; ++i) {
      kt2ij[imap[i]][jjx] = yij(p[jjx], p[imap[i]], ecm2);
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

  // store the kt2 values in the output arrays
  ev.set_y23(counter > 0 ? log10(kt2[counter - 1 - 0]) : -50.);
  ev.set_y34(counter > 1 ? log10(kt2[counter - 1 - 1]) : -50.);
  ev.set_y45(counter > 2 ? log10(kt2[counter - 1 - 2]) : -50.);
  ev.set_y56(counter > 3 ? log10(kt2[counter - 1 - 3]) : -50.);
}
