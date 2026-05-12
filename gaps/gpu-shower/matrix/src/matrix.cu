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
void calc_lome(thrust::device_vector<event>& d_events, const params& p,
               int blocks) {
  /**
   * @brief wrapper for matrix element calculation
   *
   * @param d_events device vector of events
   * @param process LHC or LEP (1 or 2)
   * @param nlo whether to calculate NLO corrections
   * @param root_s center of mass energy
   * @param asmz alpha_s at mz
   * @param blocks number of CUDA thread blocks
   * @param threads number of CUDA threads per block
   * @param pdf_name PDF set name
   */

  // number of events - can get from d_events.size()
  int n = d_events.size();

  // allocate memory for a matrix object on the device
  matrix* d_matrix;
  cudaMalloc(&d_matrix, sizeof(matrix));

  // set up the device matrix object
  debug_msg("running @matrix_setup_kernel");
  matrix_setup_kernel<<<1, 1>>>(d_matrix, p.process, p.nlo, p.root_s);
  sync_gpu_and_check("matrix_setup_kernel");

  // set up the device alpha_s calculator
  alpha_s* d_as;
  cudaMalloc(&d_as, sizeof(alpha_s));
  as_setup_kernel<<<1, 1>>>(d_as, p.asmz);
  sync_gpu_and_check("as_setup_kernel");

  // set up the pdf evaluator (for LHC processes)
  pdf_wrapper pdf(p.me2pdf);

  // LEP LO
  if ((p.process == 1) && !p.nlo) {
    lep_lo(d_events, d_matrix, blocks, p.threads);
  }

  // LEP NLO
  else if ((p.process == 1) && p.nlo) {
    lep_lo(d_events, d_matrix, blocks, p.threads);
    lep_nlo(d_events, d_matrix, d_as, blocks, p.threads);
  }

  // LHC LO - and just do NLO = LO for now
  else if ((p.process == 2) && !p.nlo) {
    lhc_lo(d_events, d_matrix, &pdf, blocks, p.threads);
  }

  else if ((p.process == 2) && p.nlo) {
    lhc_lo(d_events, d_matrix, &pdf, blocks, p.threads);
    lhc_nlo(d_events, d_matrix, d_as, &pdf, blocks, p.threads);
  }

  return;
}