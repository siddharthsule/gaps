// Attempt 1: Obtaining the Matrix Element Square on GPU

#include <curand_kernel.h>
#include <cmath>

#include "matrix.cuh"
#include "particle.cuh"
#include "vec4.cuh"

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

// Device setup function
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

__global__ void matrixSetupKernel(Matrix *matrix, double alphas, double ecms) {
  matrix->setup(alphas, ecms);
}

__global__ void loPointKernel(Matrix *matrix, Event *p_data, int n) {
  int id = blockIdx.x * blockDim.x + threadIdx.x;

  curandState state;
  curand_init(1234, id, 0, &state);

  if (id < n) {
    double ct = 2. * curand_uniform(&state) - 1.;
    double st = sqrt(1. - ct * ct);
    double phi = 2. * M_PI * curand_uniform(&state);

    int fl = curand(&state) % 5 + 1;
    double s = 91.2;
    
    // (pa - p1)^2 = -Ecms ( 1 - cos Theta), derived
    double t = -91.2 * (1 - ct);  

    double lome = matrix->ME2(fl, s, t);
    double dxs = 5. * lome * 3.89379656e8 / (8. * M_PI) / (2. * pow(91.2, 2.));

    Vec4 pa(-45.6, 0., 0., 45.6);
    Vec4 pb(-45.6, 0., 0., -45.6);
    Vec4 p1(45.6, 45.6 * st * cos(phi), 45.6 * st * sin(phi), 45.6 * ct);
    Vec4 p2(45.6, -45.6 * st * cos(phi), -45.6 * st * sin(phi), -45.6 * ct);

    Particle p[4] = {Particle(-11, pa, 0, 0), Particle(11, pb, 0, 0),
                     Particle(fl, p1, 1, 0), Particle(-fl, p2, 0, 1)};
    
    new (&p_data[id]) Event(p, dxs); 
  }
}
