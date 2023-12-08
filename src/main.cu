#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "matrix.cuh"
#include "particle.cuh"
#include "vec4.cuh"

// Jet and Event Shape Analysis
#include "durham.h"
#include "etflow.h"
#include "jetmbr.h"
#include "thrust.h"

// Dipole Shower
#include "qcd.h"
#include "shower.h"

Event *gpuLOME(const int N) {
  // Allocate memory for a Matrix object on the device
  Matrix *d_matrix;
  cudaMalloc(&d_matrix, sizeof(Matrix));

  // Set up the device Matrix object
  matrixSetupKernel<<<1, 1>>>(d_matrix, 0.118, 91.2);

  // Host and Device Variables
  Event *h_pd, *d_pd;

  h_pd = (Event *)malloc(N * sizeof(Event));
  cudaMalloc(&d_pd, N * sizeof(Event));

  // Generate the LO Matrix Elements
  loPointKernel<<<(N + 255) / 256, 256>>>(d_matrix, d_pd, N);

  // Copy the results back to the host
  cudaMemcpy(h_pd, d_pd, N * sizeof(Event), cudaMemcpyDeviceToHost);

  // Free Memory
  cudaFree(d_pd);

  return h_pd;
}

int runGenerator(const int &N, const std::string &filename = "test.yoda") {
  DAnalysis da;
  TAnalysis ta;
  MBrAnalysis ma;
  ETAnalysis ea;

  AlphaS as(91.1876, 0.1181);
  Shower sh(1., as);

  Event *events = gpuLOME(N);

  for (int i = 0; i < N; i++) {

    Event ev = events[i];

    double t = (ev.GetParton(0).GetMom() + ev.GetParton(1).GetMom()).M2();
    sh.Run(ev, t);

    //std::cout << "Event " << i << " has " << ev.GetSize() << " partons" <<
    //std::endl; for (auto& p : ev.partons) { std::cout << p.GetPid() << " ";}
    

    if (IsEventCheckValid(ev)) {
      if (ev.GetSize() > 4) {
        ta.Analyze(ev);
        da.Analyze(ev);
        ea.Analyze(ev);
        ma.Analyze(ev);
      }
    } else {
      std::cout << "Event failed validation" << std::endl;
    }

    if (i % 1000 == 0) {
      std::cout << "Running Event " << i << std::endl;
    }
  }

  ta.Finalize(filename);
  da.Finalize(filename);
  ea.Finalize(filename);
  ma.Finalize(filename);

  return 0;
}

int main(int argc, char *argv[]) {
  int N = argc > 1 ? atoi(argv[1]) : 10000;
  std::string filename = "test.yoda";
  int remover = std::remove("test.yoda");
  runGenerator(N, filename);
  return 0;
}
