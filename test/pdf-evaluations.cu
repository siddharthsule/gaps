// -----------------------------------------------------------------------------
/**
 * LHAPDF GPU PDF Evaluation for different flavours
 * ------------------------------------------------
 * This code example shows a first attempt at calculating the PDF for a bunch of
 * events, where the flavour of the event is unknown.
 *
 * The code generates random numbers for x, Q2 and flavours,
 * evaluates the PDFs and calculates the ratio between them.
 *
 * For the CPU variant this is done one by one.
 *
 * As we only have access to 2D dataset interpolation, in the GPU code, the PDFs
 * are evaluated for ALL flavours for both sets of (x, Q2) in parallel, leaving
 * us with two 2D arrays, with 11 rows (one for each flavour) and N columns, 
 * full of xf values. These arrays are then sent to a kernel, where the correct 
 * index is picked using the provided flavours and the ratio is calculated.
 */

/*
To run this code:
export INSTALL_LOC=$PWD/../gaps/lhapdf-gpu/install/
echo "INSTALL_LOC is set to: $INSTALL_LOC"
export PATH=$PATH:$INSTALL_LOC/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$INSTALL_LOC/lib
export PYTHONPATH=$PYTHONPATH:$INSTALL_LOC/lib/python3.11/site-packages
nvcc `../gaps/lhapdf-gpu/install/bin/lhapdf-config --cflags --ldflags` -o evals pdf-evaluations.cu
./evals
*/

// -----------------------------------------------------------------------------
// Includes

// cuda libraries
#include <cuda_runtime.h>
#include <curand_kernel.h>

// LHAPDF
#include "LHAPDF/CuPDF.h"
#include "LHAPDF/LHAPDF.h"

// Thrust
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

// std libraries
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>

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

void evaluate(CuPDF *pdf, int *d_fl, double *d_x, double *d_q2, double *d_xf,
              int n, int blocks, int threads) {
  /**
   * @brief Evaluate the PDFs for a given set of flavours, x and Q^2 values
   *
   * @param pdf: CuPDF pointer
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
  select_flavour<<<blocks, threads>>>(d_fl, d_flat, d_xf, n);
  cudaDeviceSynchronize();

  // Free memory
  flat_pdf_matrix.clear();

  return;
}

// -----------------------------------------------------------------------------
// Main

int main(int argc, char *argv[]) {
  // ---------------------------------------------------------------------------
  // Get N from the command line

  int N = (argc > 1) ? atoi(argv[1]) : 1000000;

  // GPU kernel settings
  int nb = 256;
  int nt = (N + nb - 1) / nb;


  // Q2 min and max
  double qmin = 1.295;
  double qmax = 100000.0;

  // ---------------------------------------------------------------------------
  // Initialize the PDF Evaluation

  LHAPDF::PDF *pdf = LHAPDF::mkPDF("CT14lo", 0);
  CuPDF *cuda_pdf = new CuPDF("CT14lo", 0);

  // ---------------------------------------------------------------------------
  // Initialize all the necessary arrays

  // Generate host arrays for x, q2 and flavours
  double *x = new double[N];
  double *q2 = new double[N];
  int *flavs = new int[N];

  // randomly generate x, q2 and flavours on the host
  for (int i = 0; i < N; ++i) {

    // x in [0, 1]
    x[i] = static_cast<double>(rand()) / RAND_MAX;

    // q2 in [qmin, qmax]
    q2[i] = qmin + (static_cast<double>(rand()) / RAND_MAX) * (qmax - qmin);

    // flavours in [1,2,3,4,5,-1,-2,-3,-4,-5,21]
    int j = rand() % 11;  // Random number between 0 and 10
    int current_flavour;
    if (j == 0) {
      current_flavour = 21;  // Gluon
    } else if (j <= 5) {
      current_flavour = j;  // Quarks
    } else {
      current_flavour = -(j - 5);  // Antiquarks
    }
    flavs[i] = current_flavour;
  }

  // Allocate memory on the device for the host arrays
  double *d_x, *d_q2;
  int *d_flavours;
  cudaMalloc(&d_x, N * sizeof(double));
  cudaMalloc(&d_q2, N * sizeof(double));
  cudaMalloc(&d_flavours, N * sizeof(int));

  // Copy the host arrays to the device
  cudaMemcpy(d_x, x, N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_q2, q2, N * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_flavours, flavs, N * sizeof(int), cudaMemcpyHostToDevice);

  // ---------------------------------------------------------------------------
  // Evaluate PDF for all flavours and calc PDF Ratio (CPU, timed)

  auto start = std::chrono::high_resolution_clock::now();

  double *xf_cpu = new double[N];
  for (int i = 0; i < N; ++i) {
    xf_cpu[i] = pdf->xfxQ2(flavs[i], x[i], q2[i]);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_cpu = end - start;

  // ---------------------------------------------------------------------------
  // Evaluate PDF for all flavours and calc PDF Ratio (GPU, timed)

  start = std::chrono::high_resolution_clock::now();

  // Allocate memory for the output arrays on the device
  double *d_xf;
  cudaMalloc(&d_xf, 11 * N * sizeof(double));  // 11 flavours

  // Evaluate PDF for all flavours
  evaluate(cuda_pdf, d_flavours, d_x, d_q2, d_xf, N, nt, nb);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_gpu = end - start;

  // ---------------------------------------------------------------------------
  // Print the time taken

  std::cout << "CPU Time: " << elapsed_cpu.count() << "s" << std::endl;
  std::cout << "GPU Time: " << elapsed_gpu.count() << "s" << std::endl;
  std::cout << "Speedup: "
            << (elapsed_cpu.count() / elapsed_gpu.count()) << "x" << std::endl;

  // ---------------------------------------------------------------------------
  // Copy the results back to the host

  double *xf = new double[N];
  cudaMemcpy(xf, d_xf, N * sizeof(double), cudaMemcpyDeviceToHost);

  // ---------------------------------------------------------------------------
  // Validation

  /**
  * THERE WON'T BE 100% MATCHES AS SOME FUNCTIONALITY IS NOT IMPLEMENTED, FOR
  * EXAMPLE SMOOTHING IN THE EDGES OF THE PDF CURVES!!!!!!!!!!!!
  */

  std::cout << "Validating results..." << std::endl;

  // Floating-point comparison with tolerance
  bool all_match = true;
  int unmatched = 0;
  for (int i = 0; i < N; ++i) {
    if (fabs(xf[i] - xf_cpu[i]) > 1e-9) {
      all_match = false;
      unmatched++;
    }
  }

  if (all_match) {
    std::cout << "All results match!" << std::endl;
  } else {
    std::cout << "Some results do not match!" << std::endl;
    std::cout << " - Matched: " << (N - unmatched) << std::endl;
    std::cout << " - Unmatched: " << unmatched << std::endl;
    std::cout << " - Percentage unmatched: "
              << (static_cast<double>(unmatched) / N) * 100.0 << "%" << std::endl;
  }

  // ---------------------------------------------------------------------------
   
  // Free memory
  delete[] x;
  delete[] q2;
  delete[] flavs;
  delete[] xf_cpu;
  delete[] xf;
  cudaFree(d_x);
  cudaFree(d_q2);
  cudaFree(d_flavours);
  cudaFree(d_xf);
  delete pdf;
  delete cuda_pdf;
  cudaDeviceReset();

  return 0;
}