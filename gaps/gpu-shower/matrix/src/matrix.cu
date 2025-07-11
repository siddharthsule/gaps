#include "matrix.cuh"

// -----------------------------------------------------------------------------
// constructor

__device__ void matrix::setup(bool nlo, double root_s) {
  this->nlo = nlo;
  this->root_s = root_s;
  this->s = pow(root_s, 2.);
  this->s_hat = s;
  this->mz2 = pow(mz, 2.);
  this->gz2 = pow(gz, 2.);
  this->alpha = 1. / 128.802;
  this->sin2tw = 0.22293;
  this->amin = 1.e-10;
  this->ye = 0.5;
  this->ze = 0.01;
  this->ws = 0.25;
}

// kernel to set up the matrix object on the device
__global__ void matrix_setup_kernel(matrix *matrix, bool nlo, double root_s) {
  matrix->setup(nlo, root_s);
}

// -----------------------------------------------------------------------------
// main

// function to generate the lo matrix elements + momenta
void calc_lome(thrust::device_vector<event> &d_events, bool nlo, double root_s,
               double asmz, int blocks, int threads) {
  /**
   * @brief wrap
   */

  // number of events - can get from d_events.size()
  int n = d_events.size();

  // allocate memory for a matrix object on the device
  matrix *d_matrix;
  cudaMalloc(&d_matrix, sizeof(matrix));

  // set up the device matrix object
  debug_msg("running @matrix_setup_kernel");
  matrix_setup_kernel<<<1, 1>>>(d_matrix, nlo, root_s);
  sync_gpu_and_check("matrix_setup_kernel");

  // set up the device alpha_s calculator
  alpha_s *d_as;
  cudaMalloc(&d_as, sizeof(alpha_s));
  as_setup_kernel<<<1, 1>>>(d_as, mz, asmz);
  sync_gpu_and_check("as_setup_kernel");

  // LEP LO
  if (!nlo) {
    debug_msg("running @lep_lo");
    lep_lo<<<blocks, threads>>>(d_matrix,
                                thrust::raw_pointer_cast(d_events.data()), n);
    sync_gpu_and_check("lep_lo");
  }

  // LEP NLO
  else {
    debug_msg("running @lep_lo and @lep_nlo");
    lep_lo<<<blocks, threads>>>(d_matrix,
                                thrust::raw_pointer_cast(d_events.data()), n);
    sync_gpu_and_check("lep_lo");
    lep_nlo<<<blocks, threads>>>(d_matrix, d_as,
                                 thrust::raw_pointer_cast(d_events.data()), n);
    sync_gpu_and_check("lep_nlo");
  }

  return;
}