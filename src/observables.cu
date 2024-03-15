#include <fstream>

#include "observables.cuh"

// -----------------------------------------------------------------------------
// Validate Events before binning

__global__ void validateEvents(Event* events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];
  ev.SetValidity(ev.Validate());

  if (!ev.GetValidity()) {
    printf("Invalid Event\n");
  }
}

// For now the event is validated before every observable, but this could be
// optimized by validating once and then passing a flag to the observables
// This is true for both C++ and CUDA, so our comparison is still valid

// -----------------------------------------------------------------------------
// Jet Rates

// Yij function Used for the Durham analysis
__device__ double Yij(const Vec4& p, const Vec4& q, double ecm2) {
  double pq = p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
  double min_pq = min(p[0], q[0]);
  double max_pq = max(pq / sqrt(p.P2() * q.P2()), -1.0);
  return 2.0 * pow(min_pq, 2) * (1.0 - min(max_pq, 1.0)) / ecm2;
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
  double kt2[maxPartons] = {-1.0};
  int counter = 0;

  // Number of partons (which will change when clustered), lower case to avoid N
  int n = ev.GetPartonSize();

  // imap will store the indices of the partons
  int imap[maxPartons];
  for (int i = 0; i < ev.GetPartonSize(); ++i) {
    imap[i] = i;
  }

  // kt2ij will store the kt2 values for each pair of partons
  double kt2ij[maxPartons][maxPartons] = {0.0};
  double dmin = 1.0;
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
    dmin = 1.0;
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
  ev.SetY23(counter > 0 ? log10(kt2[counter - 1 - 0]) : -50.0);
  ev.SetY34(counter > 1 ? log10(kt2[counter - 1 - 1]) : -50.0);
  ev.SetY45(counter > 2 ? log10(kt2[counter - 1 - 2]) : -50.0);
  ev.SetY56(counter > 3 ? log10(kt2[counter - 1 - 3]) : -50.0);
}

// -----------------------------------------------------------------------------
// Event Shapes

// Bubble sort to sort the momenta in descending order of P()
__device__ void bubbleSort(Vec4* moms, int n) {
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - i - 1; j++) {
      if (moms[j].P() < moms[j + 1].P()) {
        Vec4 temp = moms[j];
        moms[j] = moms[j + 1];
        moms[j + 1] = temp;
      }
    }
  }
}

// Thrust calculation from TASSO
__global__ void calculateThr(Event* events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];

  if (!ev.GetValidity() || ev.GetPartonSize() < 3) {
    return;
  }

  Vec4 moms[maxPartons];
  for (int i = 2; i < ev.GetSize(); ++i) {
    moms[i - 2] = ev.GetParton(i).GetMom();
  }

  bubbleSort(moms, maxPartons);

  double momsum = 0.0;
  for (int i = 0; i < ev.GetPartonSize(); ++i) {
    momsum += moms[i].P();
  }

  double thr = 0.0;
  Vec4 t_axis = Vec4();

  for (int k = 1; k < ev.GetPartonSize(); ++k) {
    for (int j = 0; j < k; ++j) {
      Vec4 tmp_axis = moms[j].Cross(moms[k]);
      Vec4 p_thrust = Vec4();
      Vec4 p_combin[4];

      for (int i = 0; i < ev.GetPartonSize(); ++i) {
        if (i != j && i != k) {
          if (moms[i].Dot(tmp_axis) >= 0) {
            p_thrust = p_thrust + moms[i];
          } else {
            p_thrust = p_thrust - moms[i];
          }
        }
      }

      p_combin[0] = (p_thrust + moms[j] + moms[k]);
      p_combin[1] = (p_thrust + moms[j] - moms[k]);
      p_combin[2] = (p_thrust - moms[j] + moms[k]);
      p_combin[3] = (p_thrust - moms[j] - moms[k]);

      for (int i = 0; i < 4; ++i) {
        double temp = p_combin[i].P();
        if (temp > thr) {
          thr = temp;
          t_axis = p_combin[i];
        }
      }
    }
  }

  thr /= momsum;
  thr = 1.0 - thr;

  t_axis = t_axis / (t_axis).P();
  if (t_axis[3] < 0) {
    t_axis = t_axis * -1.0;
  }

  if (thr < 1e-12) {
    thr = -5.0;
  }

  ev.SetThr(thr);
  ev.SetTAxis(t_axis);
}

// Jet Mass and Broadening
__global__ void CalculateJetMBr(Event* events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];

  if (!ev.GetValidity() || ev.GetPartonSize() < 3) {
    return;
  }

  Vec4 moms[maxPartons];
  for (int i = 2; i < ev.GetSize(); ++i) {
    moms[i - 2] = ev.GetParton(i).GetMom();
  }

  double momsum = 0.0;
  for (int i = 0; i < ev.GetSize(); ++i) {
    momsum += moms[i].P();
  }

  Vec4 p_with, p_against;
  int n_with = 0, n_against = 0;
  double e_vis = 0.0, broad_with = 0.0, broad_against = 0.0,
         broad_denominator = 0.0;

  for (int i = 0; i < ev.GetPartonSize(); ++i) {
    double mo_para = moms[i].Dot(ev.GetTAxis());
    double mo_perp = (moms[i] - (ev.GetTAxis() * mo_para)).P();
    double enrg = moms[i].P();

    e_vis += enrg;
    broad_denominator += 2.0 * enrg;

    if (mo_para > 0.0) {
      p_with = p_with + moms[i];
      broad_with += mo_perp;
      n_with++;
    } else if (mo_para < 0.0) {
      p_against = p_against + moms[i];
      broad_against += mo_perp;
      n_against++;
    } else {
      p_with = p_with + (moms[i] * 0.5);
      p_against = p_against + (moms[i] * 0.5);
      broad_with += 0.5 * mo_perp;
      broad_against += 0.5 * mo_perp;
      n_with++;
      n_against++;
    }
  }

  double e2_vis = e_vis * e_vis;

  double mass2_with = fabs(p_with.M2() / e2_vis);
  double mass2_against = fabs(p_against.M2() / e2_vis);

  double mass_with = sqrt(mass2_with);
  double mass_against = sqrt(mass2_against);

  broad_with /= broad_denominator;
  broad_against /= broad_denominator;

  double mH = fmax(mass_with, mass_against);
  double mL = fmin(mass_with, mass_against);

  double bW = fmax(broad_with, broad_against);
  double bN = fmin(broad_with, broad_against);

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
// Analysis

// Fill theHistograms (Atomically!)
__global__ void fillHistos(Analysis* an, Event* events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= N) {
    return;
  }

  Event& ev = events[idx];

  an->hists[0].Fill(ev.GetY23(), ev.GetDxs());
  an->hists[1].Fill(ev.GetY34(), ev.GetDxs());
  an->hists[2].Fill(ev.GetY45(), ev.GetDxs());
  an->hists[3].Fill(ev.GetY56(), ev.GetDxs());
  an->hists[4].Fill(ev.GetThr(), ev.GetDxs());
  an->hists[5].Fill(ev.GetThr(), ev.GetDxs());
  an->hists[6].Fill(ev.GetHJM(), ev.GetDxs());
  an->hists[7].Fill(ev.GetLJM(), ev.GetDxs());
  an->hists[8].Fill(ev.GetWJB(), ev.GetDxs());
  an->hists[9].Fill(ev.GetNJB(), ev.GetDxs());

  atomicAdd(&an->wtot, ev.GetDxs());
  atomicAdd(&an->ntot, 1.0);
}

// Run the above kernels
void doAnalysis(thrust::device_vector<Event>& d_events, std::string filename) {
  /**
   * Only Place for a Host Object
   * ----------------------------
   *
   * While most of the work is done on the device, one cannot directly write
   * to a file from the device. Therefore, we will create a host object to
   * store the histograms, and then copy the results back to the host for
   * writing to file.
   */

  // Device Analysis Object
  Analysis *h_an, *d_an;

  // Allocate memory for the device analysis object
  h_an = new Analysis();
  cudaMalloc(&d_an, sizeof(Analysis));
  cudaMemcpy(d_an, h_an, sizeof(Analysis), cudaMemcpyHostToDevice);

  // Get Event Data
  int N = d_events.size();
  Event* d_events_ptr = thrust::raw_pointer_cast(d_events.data());

  // Validate the Events
  validateEvents<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("validateEvents");

  // Calculare the Observables
  doCluster<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("doCluster");

  calculateThr<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("calculateThr");

  CalculateJetMBr<<<(N + 255) / 256, 256>>>(d_events_ptr, N);
  syncGPUAndCheck("CalculateJetMBr");

  // Do the Analysis
  fillHistos<<<(N + 255) / 256, 256>>>(d_an, d_events_ptr, N);
  syncGPUAndCheck("fillHistos");

  // Copy the results back to the host
  cudaMemcpy(h_an, d_an, sizeof(Analysis), cudaMemcpyDeviceToHost);

  // Normalize the histograms
  for (auto& hist : h_an->hists) {
    hist.ScaleW(1.0 / h_an->ntot);
  }

  // Remove existing file
  std::remove(filename.c_str());

  // Write the histograms to file
  Write(h_an->hists[0], "/gaps/log10y23\n", filename);
  Write(h_an->hists[1], "/gaps/log10y34\n", filename);
  Write(h_an->hists[2], "/gaps/log10y45\n", filename);
  Write(h_an->hists[3], "/gaps/log10y56\n", filename);
  Write(h_an->hists[4], "/gaps/tvalue\n", filename);
  Write(h_an->hists[5], "/gaps/tzoomd\n", filename);
  Write(h_an->hists[6], "/gaps/hjm\n", filename);
  Write(h_an->hists[7], "/gaps/ljm\n", filename);
  Write(h_an->hists[8], "/gaps/wjb\n", filename);
  Write(h_an->hists[9], "/gaps/njb\n", filename);

  // Clean up
  delete h_an;
  cudaFree(d_an);
}