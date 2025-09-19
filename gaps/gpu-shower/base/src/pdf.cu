#include "pdf.cuh"

// -----------------------------------------------------------------------------
// After evaluating all flavours (below), filter data to retrieve right value

__global__ void select_flavour(int *fl, double *xf_in, double *xf_out, int n) {
  /**
   * @brief Filter the data to retrieve the right value
   *
   * @param fl: Flavour array
   * @param xf_in: Input array
   * @param xf_out: Output array
   * @param n: Number of elements
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------

  // Convert Flavour to Row
  int fl_row = (fl[idx] == 21) ? 0
               : (fl[idx] > 0) ? fl[idx]
               : (fl[idx] < 0) ? (-fl[idx] + 5)
                               : -1;

  // Handle invalid flavours
  if (fl_row == -1) {
    xf_out[idx] = 0.0;
    return;
  }

  // Calculate the right index for the flat array
  // Formula: [row][col] -> [row * n + col]
  int idx_flat = fl_row * n + idx;

  // Get the right value
  xf_out[idx] = xf_in[idx_flat];

  return;
}

// -----------------------------------------------------------------------------
// wrapper function that calls the kernel

void pdf_wrapper::evaluate(int *d_fl, double *d_x, double *d_q2, double *d_xf,
                           int n, int blocks, int threads) {
  /**
   * @brief Evaluate the PDFs for a given set of flavours, x and Q^2 values
   *
   * @param d_fl: Flavour array
   * @param d_x: x values
   * @param d_q2: Q^2 values
   * @param d_xf: Output array
   * @param n: Number of elements
   */

  // Allocate a flat 2D array on the device
  // Formula: [row][col] -> [row * n + col]
  thrust::device_vector<double> flat_pdf_matrix(11 * n, 0.);
  double *d_flat = thrust::raw_pointer_cast(flat_pdf_matrix.data());

  // evaluate PDF for all quarks, antiquarks and gluon
  for (int i = 0; i < 11; ++i) {
    // Calculate the correct offset for each row
    double *current_row = d_flat + i * n;

    // Get the flavour - 0 is for gluon, 1-5 for quarks, 6-10 for antiquarks
    int current_flavour;
    if (i == 0) {
      current_flavour = 21;  // Gluon
    } else if (i <= 5) {
      current_flavour = i;  // Quarks
    } else {
      current_flavour = -(i - 5);  // Antiquarks
    }

    ///////////////////////////////////////////////////////
    // Evaluate PDF for a flavour - EXTERNAL
    pdf->xfxQ2(d_x, d_q2, current_row, current_flavour, n);
    ///////////////////////////////////////////////////////
  }

  // After evaluating all flavours, filter data to retrieve right value
  debug_msg("running @select_flavour");
  select_flavour<<<blocks, threads>>>(d_fl, d_flat, d_xf, n);
  sync_gpu_and_check("select_flavour");

  // Free memory
  flat_pdf_matrix.clear();

  return;
}