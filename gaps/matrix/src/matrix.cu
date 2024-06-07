#include "matrix.cuh"
#include "prng.cuh"

// host constructor
matrix::matrix(double alphas, double ecms)
    : alphas(alphas),
      ecms(ecms),
      mz2(pow(91.1876, 2.)),
      gz2(pow(2.4952, 2.)),
      alpha(1. / 128.802),
      sin2tw(0.22293),
      amin(1.e-10),
      ye(0.5),
      ze(0.01),
      ws(0.25) {}

// device setup function - default values in matrix.cuh
__device__ void matrix::setup(double alphas, double ecms) {
  this->alphas = alphas;
  this->ecms = ecms;
  this->mz2 = pow(91.1876, 2.);
  this->gz2 = pow(2.4952, 2.);
  this->alpha = 1. / 128.802;
  this->sin2tw = 0.22293;
  this->amin = 1.e-10;
  this->ye = 0.5;
  this->ze = 0.01;
  this->ws = 0.25;
}

// me^2 formula
__device__ double matrix::me2(int fl, double s, double t) {
  double qe = -1.;
  double ae = -0.5;
  double ve = ae - 2. * qe * sin2tw;
  double qf = (fl == 2 || fl == 4) ? 2. / 3. : -1. / 3.;
  double af = (fl == 2 || fl == 4) ? 0.5 : -0.5;
  double vf = af - 2. * qf * sin2tw;
  double kappa = 1. / (4. * sin2tw * (1. - sin2tw));
  double chi1 = kappa * s * (s - mz2) / (pow(s - mz2, 2.) + gz2 * mz2);
  double chi2 = pow(kappa * s, 2.) / (pow(s - mz2, 2.) + gz2 * mz2);
  double term1 = (1. + pow(1. + 2. * t / s, 2.)) *
                 (pow(qf * qe, 2.) + 2. * (qf * qe * vf * ve) * chi1 +
                  (ae * ae + ve * ve) * (af * af + vf * vf) * chi2);
  double term2 = (1. + 2. * t / s) * (4. * qe * qf * ae * af * chi1 +
                                      8. * ae * ve * af * vf * chi2);
  return pow(4. * M_PI * alpha, 2.) * 3. * (term1 + term2);
}

// kernel to set up the matrix object on the device
__global__ void matrix_setup_kernel(matrix *matrix, double e) {
  matrix->setup(asmz, e);
}

// kernel to generate the event
__global__ void lo_point_kernel(matrix *matrix, event *events, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event &ev = events[idx];

  // Set seed and random number
  unsigned long seed = ev.get_seed();
  double rand = ev.get_rand();

  // First Time Only - Generate the random number
  update_rng(seed, rand);
  ev.set_seed(seed);
  ev.set_rand(rand);
  if (idx == 0) {
    printf("Seed: %lu, Rand: %f, First Gen \n", seed, rand);
  }

  int fl = static_cast<int>(rand) % 5 + 1;
  update_rng(seed, rand);
  ev.set_seed(seed);
  ev.set_rand(rand);
  if (idx == 0) {
    printf("Seed: %lu, Rand: %f, Flavour \n", seed, rand);
  }

  double ct = 2. * rand - 1.;
  update_rng(seed, rand);
  ev.set_seed(seed);
  ev.set_rand(rand);
  if (idx == 0) {
    printf("Seed: %lu, Rand: %f, Cos Theta \n", seed, rand);
  }
  double st = sqrt(1. - ct * ct);
  double phi = 2. * M_PI * rand;
  update_rng(seed, rand);
  ev.set_seed(seed);
  ev.set_rand(rand);
  if (idx == 0) {
    printf("Seed: %lu, Rand: %f, Phi \n", seed, rand);
  }

  double p0 = matrix->get_ecms() / 2.;  // need to use get because outside class

  vec4 pa(p0, 0., 0., p0);
  vec4 pb(p0, 0., 0., -p0);
  vec4 p1(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);
  vec4 p2(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);

  double lome = matrix->me2(fl, (pa + pb).m2(), (pa - p1).m2());

  // calculate the differential cross section
  // 5 = 5 flavours (?)
  // 3.89379656e8 = convert from ge_v^-2 to pb
  // 8 pi = standard phase space factor
  // pow(matrix->get_ecms(), 2.) = center of mass energy squared, s
  double dxs = 5. * lome * 3.89379656e8 / (8. * M_PI) /
               (2. * pow(matrix->get_ecms(), 2.));

  parton p[4] = {parton(-11, -pa, 0, 0), parton(11, -pb, 0, 0),
                 parton(fl, p1, 1, 0), parton(-fl, p2, 0, 1)};

  // set the partons
  for (int i = 0; i < 4; i++) {
    ev.set_parton(i, p[i]);
  }

  // set the me params
  ev.set_dxs(dxs);
  ev.set_hard(4);

  // set the random seed
  ev.set_seed(seed);
  ev.set_rand(rand);
}

// function to generate the lo matrix elements + momenta
void calc_lome(thrust::device_vector<event> &d_events, double e) {
  // number of events - can get from d_events.size()
  int n = d_events.size();

  // allocate memory for a matrix object on the device
  matrix *d_matrix;
  cudaMalloc(&d_matrix, sizeof(matrix));

  // set up the device matrix object
  debug_msg("running @matrix_setup_kernel");
  matrix_setup_kernel<<<1, 1>>>(d_matrix, e);
  sync_gpu_and_check("matrix_setup_kernel");

  // generate the lo matrix elements
  debug_msg("running @lo_point_kernel");
  lo_point_kernel<<<(n + 255) / 256, 256>>>(
      d_matrix, thrust::raw_pointer_cast(d_events.data()), n);
  sync_gpu_and_check("lo_point_kernel");

  // free memory
  cudaFree(d_matrix);

  return;
}
