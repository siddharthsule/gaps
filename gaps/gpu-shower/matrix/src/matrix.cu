#include "matrix.cuh"

// -----------------------------------------------------------------------------
// constructor

__device__ void matrix::setup(int process, bool nlo, double root_s) {
  this->process = process;
  this->nlo = nlo;
  this->root_s = root_s;
  this->ws = 0.25;
}

// kernel to set up the matrix object on the device
__global__ void matrix_setup_kernel(matrix* matrix, int process, bool nlo,
                                    double root_s) {
  matrix->setup(process, nlo, root_s);
}

// -----------------------------------------------------------------------------
// main

// function to generate the lo matrix elements + momenta
void calc_lome(thrust::device_vector<event>& d_events, int process, bool nlo,
               double root_s, double asmz, int blocks, int threads) {
  /**
   * @brief wrap
   */

  // number of events - can get from d_events.size()
  int n = d_events.size();

  // allocate memory for a matrix object on the device
  matrix* d_matrix;
  cudaMalloc(&d_matrix, sizeof(matrix));

  // set up the device matrix object
  debug_msg("running @matrix_setup_kernel");
  matrix_setup_kernel<<<1, 1>>>(d_matrix, process, nlo, root_s);
  sync_gpu_and_check("matrix_setup_kernel");

  // set up the device alpha_s calculator
  alpha_s* d_as;
  cudaMalloc(&d_as, sizeof(alpha_s));
  as_setup_kernel<<<1, 1>>>(d_as, asmz);
  sync_gpu_and_check("as_setup_kernel");

  // LEP LO
  if ((process == 1) && !nlo) {
    lep_lo(d_events, d_matrix, blocks, threads);
  }

  // LEP NLO
  else if ((process == 1) && nlo) {
    lep_lo(d_events, d_matrix, blocks, threads);
    lep_nlo(d_events, d_matrix, d_as, blocks, threads);
  }

  // LHC LO - and just do NLO = LO for now
  else if ((process == 2) && !nlo) {
    lhc_lo(d_events, d_matrix, blocks, threads);
  }

  else if ((process == 2) && nlo) {
    lhc_lo(d_events, d_matrix, blocks, threads);
    lhc_nlo(d_events, d_matrix, d_as, blocks, threads);
  }

  return;
}