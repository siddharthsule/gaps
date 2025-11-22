#include "matrix.cuh"

// -----------------------------------------------------------------------------
// Utility Kernels

__global__ void multiply_dxs_by_pdf_lo(event* events, int n, double* xf) {
  /**
   * @brief Multiply the cross section of each event by the PDFs
   *
   * @param events array of event records
   * @param n number of events
   * @param xf_a array of PDF values
   *
   * To improve in a future version:
   * - This is a copy of the kernel in lhc_nlo.cu, with a slightly different
   *   name, later to do all different calculations in device functions
   *   and use one kernel
   *
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Matrix Preamble
  event& ev = events[idx];
  // ---------------------------------------------

  // Multiply by PDF
  double pdf = xf[idx];
  ev.set_dxs(ev.get_dxs() * pdf);
}

// -----------------------------------------------------------------------------

__global__ void lhc_lo_no_pdf(event* events, int n, matrix* matrix, int* fl_a,
                              int* fl_b, double* x_a, double* x_b,
                              double* mu2) {
  /**
   * @brief generate the leading order interaction p p -> e+ e-
   *
   * The Leading order differential cross section for p p -> e+ e- is
   * given by:
   *
   * dsigma = drho1 drho2 drho3 drho4 * etaf(eta_a, s_hat) etaf(eta_b, s_hat)
   *          * (1/s_hat) * (I+ - I-)((s_hat - mz2)2 + mz2 gz2)/(mz gz)
   *          * ln(s/s_hat) * 1/(2 s_hat) 1/(8pi) * 2 * 1/N_c
   *          * 1/4 Sum_spins |M|^2
   *
   * @param ev the event object
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  event& ev = events[idx];
  // ---------------------------------------------

  // LHC: Bias Flavour based on contribution (based on cross section)
  // Flavours:    [down,       up,         strange,    charm,      bottom    ]
  // Percentages: [0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279]
  // Cumulatives: [0.36677679, 0.80186718, 0.91861688, 0.97129721, 1.00000000]
  double pd[5] = {0.36677679, 0.43509039, 0.11674970, 0.05268033, 0.02870279};
  double pc[5] = {0.36677679, 0.80186718, 0.91861688, 0.97129721, 1.00000000};

  // For unbiased Flavour Production
  // double pd[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
  // double pc[5] = {0.2, 0.4, 0.6, 0.8, 1.0};

  // Generate Flavor
  int fl;
  double r = ev.gen_random();
  for (int i = 0; i < 5; ++i) {
    if (r < pc[i]) {
      fl = i + 1;  // Flavors are 1-based
      break;
    }
  }

  // Four random numbers for s_hat, y, cos(theta) and phi
  double rho_1 = ev.gen_random();
  double rho_2 = ev.gen_random();
  double rho_3 = ev.gen_random();
  double rho_4 = ev.gen_random();

  // Generate s_hat
  double I_a = atan(((mz_cut_a * mz_cut_a) - mz * mz) / (mz * gz));
  double I_b = atan(((mz_cut_b * mz_cut_b) - mz * mz) / (mz * gz));
  double I = I_a + (I_b - I_a) * rho_1;
  double s_hat = (mz * gz * tan(I)) + (mz * mz);

  // generate y
  double limit =
      fmin(100., 0.5 * log((matrix->root_s * matrix->root_s) / s_hat));
  double y = -limit + 2. * limit * rho_2;

  // Generate cos(theta), phi
  double ct = -1. + 2. * rho_3;
  double phi = 2. * M_PI * rho_4;

  // Kinematics
  // Momentum Fractions
  double eta_a = sqrt(s_hat) / matrix->root_s * exp(y);
  double eta_b = sqrt(s_hat) / matrix->root_s * exp(-y);
  // Quark / AntiQuark Momentum
  vec4 pa =
      vec4(eta_a * 0.5 * matrix->root_s, 0., 0., eta_a * 0.5 * matrix->root_s);
  vec4 pb =
      vec4(eta_b * 0.5 * matrix->root_s, 0., 0., -eta_b * 0.5 * matrix->root_s);
  // Z Momentum
  vec4 q = pa + pb;
  // Z Rest Frame
  vec4 q_rest = q.boost(q);
  // Electron Positron in Rest Frame
  double p0 = 0.5 * q_rest.m();
  double st = sqrt(1. - ct * ct);  // sin(theta)
  vec4 pr1 = vec4(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);
  vec4 pr2 = vec4(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);
  // Boost back to Lab Frame
  vec4 p1 = q.boost_back(pr1);
  vec4 p2 = q.boost_back(pr2);

  // Two Possibilities: P(q) P(qbar) <-> P(qbar) P(q)
  // In this case, swap momentum and momentum fraction
  if (ev.gen_random() < 0.5) {
    vec4 temp_p = pa;
    pa = pb;
    pb = temp_p;
    double temp_eta = eta_a;
    eta_a = eta_b;
    eta_b = temp_eta;
  }

  // Generate the particles
  particle p[4] = {particle(), particle(), particle(), particle()};
  p[0] = particle(fl, pa, 0, 1, eta_a);
  p[1] = particle(-fl, pb, 1, 0, eta_b);
  p[2] = particle(-11, p1, 0, 0);
  p[3] = particle(11, p2, 0, 0);

  // Calculate the Matrix Element
  double lome = matrix->me2(fl, (pa + pb).m2(), (pa - p1).m2());
  lome *= 1 / k_nc;  // Three possible initial colour states

  // Store PDF Calculation details
  fl_a[idx] = fl;
  fl_b[idx] = -fl;
  x_a[idx] = eta_a;
  x_b[idx] = eta_b;
  mu2[idx] = s_hat;

  // Calculate the Cross Section
  double dxs;
  dxs = (1. / s_hat);
  dxs *= (I_b - I_a);
  dxs *= (pow((s_hat - mz * mz), 2) + mz * mz * gz * gz) / (mz * gz);
  dxs *= log((matrix->root_s * matrix->root_s) / s_hat);
  dxs *= (1. / (2. * s_hat)) * (1. / (8. * M_PI)) * lome;
  dxs *= GeV_minus_2_to_pb;  // 5 flavours + units
  dxs /= pd[abs(fl) - 1];    // Flavour Selection
  dxs *= 2.;                 // Two Possible Orientations

  // Set the particles
  for (int i = 0; i < 4; i++) {
    ev.set_particle(i, p[i]);
  }
  ev.set_hard(4);

  // Store the Matrix Element and Cross Section
  ev.set_me2(lome);
  ev.set_dxs(dxs);
}

// -----------------------------------------------------------------------------
// Wrapper to do the Kernels and External LHAPDF calculation

void lhc_lo(thrust::device_vector<event>& dv_events, matrix* matrix, int blocks,
            int threads) {
  /**
   * @brief run the hadronic_dxs conversion
   *
   * @param dv_events device vector of event records
   * @param proces LHC, DIS or LEP
   * @param root_s partonic centre of mass energy
   * @param blocks number of blocks to use
   * @param threads number of threads per block
   *
   */

  // set up the pdf evaluator
  pdf_wrapper pdf;

  // use a pointer to the device events
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n = dv_events.size();

  int* d_fl_a;
  cudaMalloc(&d_fl_a, n * sizeof(int));
  int* d_fl_b;
  cudaMalloc(&d_fl_b, n * sizeof(int));
  double* d_x_a;
  cudaMalloc(&d_x_a, n * sizeof(double));
  double* d_x_b;
  cudaMalloc(&d_x_b, n * sizeof(double));
  double* d_mu2;
  cudaMalloc(&d_mu2, n * sizeof(double));
  double* d_xf_a;
  cudaMalloc(&d_xf_a, n * sizeof(double));
  double* d_xf_b;
  cudaMalloc(&d_xf_b, n * sizeof(double));

  debug_msg("running @lhc_lo_no_pdf");
  lhc_lo_no_pdf<<<blocks, threads>>>(d_events, n, matrix, d_fl_a, d_fl_b, d_x_a,
                                     d_x_b, d_mu2);
  sync_gpu_and_check("lhc_lo_no_pdf");

  // LHC
  pdf.evaluate(d_fl_a, d_x_a, d_mu2, d_xf_a, n, blocks, threads);
  pdf.evaluate(d_fl_b, d_x_b, d_mu2, d_xf_b, n, blocks, threads);
  multiply_dxs_by_pdf_lo<<<blocks, threads>>>(d_events, n, d_xf_a);
  multiply_dxs_by_pdf_lo<<<blocks, threads>>>(d_events, n, d_xf_b);

  // free memory
  cudaFree(d_x_a);
  cudaFree(d_x_b);
  cudaFree(d_fl_a);
  cudaFree(d_fl_b);
  cudaFree(d_mu2);
  cudaFree(d_xf_a);
  cudaFree(d_xf_b);
  // pdf_wrapper destructs when function ends

  return;
}
