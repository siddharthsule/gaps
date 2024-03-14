#include "matrix.cuh"

// Host constructor
Matrix::Matrix(double alphas, double ecms)
    : alphas(alphas),
      ecms(ecms),
      MZ2(pow(91.1876, 2.)),
      GZ2(pow(2.4952, 2.)),
      alpha(1. / 128.802),
      sin2tw(0.22293),
      amin(1.e-10),
      ye(0.5),
      ze(0.01),
      ws(0.25) {}

// Device setup function - Default Values in matrix.cuh
__device__ void Matrix::setup(double alphas, double ecms) {
  this->alphas = alphas;
  this->ecms = ecms;
  this->MZ2 = pow(91.1876, 2.);
  this->GZ2 = pow(2.4952, 2.);
  this->alpha = 1. / 128.802;
  this->sin2tw = 0.22293;
  this->amin = 1.e-10;
  this->ye = 0.5;
  this->ze = 0.01;
  this->ws = 0.25;
}

// ME^2 Formula
__device__ double Matrix::ME2(int fl, double s, double t) {
  double qe = -1.;
  double ae = -0.5;
  double ve = ae - 2. * qe * sin2tw;
  double qf = (fl == 2 || fl == 4) ? 2. / 3. : -1. / 3.;
  double af = (fl == 2 || fl == 4) ? 0.5 : -0.5;
  double vf = af - 2. * qf * sin2tw;
  double kappa = 1. / (4. * sin2tw * (1. - sin2tw));
  double chi1 = kappa * s * (s - MZ2) / (pow(s - MZ2, 2.) + GZ2 * MZ2);
  double chi2 = pow(kappa * s, 2.) / (pow(s - MZ2, 2.) + GZ2 * MZ2);
  double term1 = (1. + pow(1. + 2. * t / s, 2.)) *
                 (pow(qf * qe, 2.) + 2. * (qf * qe * vf * ve) * chi1 +
                  (ae * ae + ve * ve) * (af * af + vf * vf) * chi2);
  double term2 = (1. + 2. * t / s) * (4. * qe * qf * ae * af * chi1 +
                                      8. * ae * ve * af * vf * chi2);
  return pow(4. * M_PI * alpha, 2.) * 3. * (term1 + term2);
}

// Kernel to set up the Matrix object on the device
__global__ void matrixSetupKernel(Matrix *matrix, double E) { matrix->setup(asmz, E); }

// Kernel to generate the Event
__global__ void loPointKernel(Matrix *matrix, Event *events, int N) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  curandState state;
  curand_init(clock64(), idx, 0, &state);

  if (idx >= N) {
    return;
  }

  double ct = 2. * curand_uniform(&state) - 1.;
  double st = sqrt(1. - ct * ct);
  double phi = 2. * M_PI * curand_uniform(&state);

  int fl = curand(&state) % 5 + 1;
  double p0 = matrix->GetECMS() / 2.;  // Need to use Get because outside class

  Vec4 pa(p0, 0., 0., p0);
  Vec4 pb(p0, 0., 0., -p0);
  Vec4 p1(p0, p0 * st * cos(phi), p0 * st * sin(phi), p0 * ct);
  Vec4 p2(p0, -p0 * st * cos(phi), -p0 * st * sin(phi), -p0 * ct);

  double lome = matrix->ME2(fl, (pa + pb).M2(), (pa - p1).M2());
  double dxs = 5. * lome * 3.89379656e8 / (8. * M_PI) /
               (2. * pow(matrix->GetECMS(), 2.));

  Parton p[4] = {Parton(-11, -pa, 0, 0), Parton(11, -pb, 0, 0),
                 Parton(fl, p1, 1, 0), Parton(-fl, p2, 0, 1)};

  Event &ev = events[idx];

  // Set the Partons
  for (int i = 0; i < 4; i++) {
    ev.SetParton(i, p[i]);
  }

  // Set the ME Params
  ev.SetDxs(dxs);
  ev.SetHard(4);
}

// Function to generate the LO Matrix Elements + Momenta
void calcLOME(thrust::device_vector<Event> &d_events, double E) {
  // Number of Events - Can get from d_events.size()
  int N = d_events.size();

  // Allocate memory for a Matrix object on the device
  Matrix *d_matrix;
  cudaMalloc(&d_matrix, sizeof(Matrix));

  // Set up the device Matrix object
  DEBUG_MSG("Running @matrixSetupKernel");
  matrixSetupKernel<<<1, 1>>>(d_matrix, E);
  syncGPUAndCheck("matrixSetupKernel");

  // Generate the LO Matrix Elements
  DEBUG_MSG("Running @loPointKernel");
  loPointKernel<<<(N + 255) / 256, 256>>>(
      d_matrix, thrust::raw_pointer_cast(d_events.data()), N);
  syncGPUAndCheck("loPointKernel");

  // Free Memory
  cudaFree(d_matrix);

  return;
}
