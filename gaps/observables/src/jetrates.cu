#include "jetrates.cuh"

// jet rates

// yij function used for the durham analysis
__device__ double yij(const vec4& p, const vec4& q, double ecm2) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = min(p[0], q[0]);
  double max_pq = max(pq / sqrt(p.p2() * q.p2()), -1.);
  return 2. * pow(min_pq, 2) * (1. - min(max_pq, 1.)) / ecm2;
}

// durham clustering algorithm
__global__ void do_cluster(event* events, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event& ev = events[idx];

  if (!ev.get_validity()) {
    return;
  }

  /**
   * on the size of arrays during the clustering process:
   *
   * the number of partons in the event is not known at compile time, so we
   * cannot use a fixed size array. we could use a dynamic array, but that
   * would require a lot of memory management, and we would have to use
   * malloc and free. instead, we will use a fixed size array, and we will
   * assume that the number of partons will not exceed max_partons. this is
   * not a great solution, but ok for now.
   */

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

  // number of partons (which will change when clustered), to avoid n
  int n_partons = ev.get_parton_size();

  // imap will store the indices of the partons
  int imap[max_partons];
  for (int i = 0; i < ev.get_parton_size(); ++i) {
    imap[i] = i;
  }

  // kt2ij will store the kt2 values for each pair of partons
  double kt2ij[max_partons][max_partons] = {0.};
  double dmin = 1.;
  int ii = 0, jj = 0;
  for (int i = 0; i < n_partons; ++i) {
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
  while (n_partons > 2) {
    --n_partons;
    kt2[counter] = dmin;
    counter++;
    int jjx = imap[jj];
    p[jjx] = p[jjx] + p[imap[ii]];
    for (int i = ii; i < n_partons; ++i) {
      imap[i] = imap[i + 1];
    }
    for (int j = 0; j < jj; ++j) {
      kt2ij[jjx][imap[j]] = yij(p[jjx], p[imap[j]], ecm2);
    }
    for (int i = jj + 1; i < n_partons; ++i) {
      kt2ij[imap[i]][jjx] = yij(p[jjx], p[imap[i]], ecm2);
    }
    dmin = 1.;
    for (int i = 0; i < n_partons; ++i) {
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