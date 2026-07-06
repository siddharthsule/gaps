#include "matrix.cuh"

// -----------------------------------------------------------------------------
// constructor

__device__ void matrix::setup(int process, bool nlo, double root_s) {
  /**
   * @brief construct the matrix element generator
   *
   * stores the run configuration used by all matrix element functions: which
   * hard process to generate (lep or lhc), whether to run at lo or nlo, and
   * the collider center-of-mass energy.
   *
   * @param process the hard process: 1 = lep (e+e- -> qqbar), 2 = lhc (pp -> Z)
   * @param nlo whether to generate at nlo (true) or lo (false)
   * @param root_s the collider center-of-mass energy
   */
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
void run_matrix(thrust::device_vector<event>& d_events, const params& p,
                int blocks) {
  /**
   * @brief wrapper for matrix element calculation
   *
   * @param d_events device vector of events
   * @param p run parameters
   * @param blocks number of CUDA thread blocks
   */

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
  std::string pdf_name_corrected = p.me2pdf != "Null"
                                       ? p.me2pdf
                                       : p.nlo ? "NNPDF40MC_nlo_as_01180"
                                               : "NNPDF40MC_lo_as_01180";
  pdf_wrapper pdf(pdf_name_corrected);

  // LEP LO
  if ((p.process == 1) && !p.nlo) {
    lep_lo(d_events, d_matrix, blocks, p.threads);
  }

  // LEP NLO
  else if ((p.process == 1) && p.nlo) {
    lep_lo(d_events, d_matrix, blocks, p.threads);
    lep_nlo(d_events, d_matrix, d_as, blocks, p.threads);
  }

  // LHC LO
  else if ((p.process == 2) && !p.nlo) {
    lhc_lo(d_events, d_matrix, &pdf, blocks, p.threads);
  }

  // LHC NLO
  else if ((p.process == 2) && p.nlo) {
    lhc_lo(d_events, d_matrix, &pdf, blocks, p.threads);
    lhc_nlo(d_events, d_matrix, d_as, &pdf, blocks, p.threads);
  }

  cudaFree(d_matrix);
  cudaFree(d_as);
}