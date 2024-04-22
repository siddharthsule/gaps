#include "jetrates.cuh"

// Jet Rates

// Yij function Used for the Durham analysis
__device__ double Yij(const Vec4& p, const Vec4& q, double ecm2) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = min(p[0], q[0]);
  double max_pq = max(pq / sqrt(p.P2() * q.P2()), -1.);
  return 2. * pow(min_pq, 2) * (1. - min(max_pq, 1.)) / ecm2;
}

// Durham Clustering Algorithm
__global__ void doCluster(Event* events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];

  if (!ev.GetValidity()) {
    return;
  }

  /**
   * On the size of arrays during the clustering process:
   *
   * The number of partons in the event is not known at compile time, so we
   * cannot use a fixed size array. We could use a dynamic array, but that
   * would require a lot of memory management, and we would have to use
   * malloc and free. Instead, we will use a fixed size array, and we will
   * assume that the number of partons will not exceed maxPartons. This is
   * not a great solution, but ok for now.
   */

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