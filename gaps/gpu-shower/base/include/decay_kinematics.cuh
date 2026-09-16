#ifndef decay_kinematics_cuh_
#define decay_kinematics_cuh_

#include "base.cuh"
#include "database.cuh"
#include "event.cuh"
#include "vec4.cuh"

// Decay kinematics, shared by cluster fission, gluon splitting, cluster decay
// and hadron decay: the Kallen function, the rest frame momentum it gives, and
// the throws that turn a parent momentum and a set of child masses into lab
// frame momenta.

// Max trials for hit-or-miss
const int max_phase_space_trials = 10000;

// The Kallen function of three squared masses
__device__ inline double kallen_lambda(double a, double b, double c) {
  return a * a + b * b + c * c - 2. * a * b - 2. * a * c - 2. * b * c;
}

// The momentum each product carries in the rest frame of a two-body decay
// m0 -> m1 m2, which supplies the two-body phase space
__device__ inline double p_star(double m0, double m1, double m2) {
  return 0.5 * sqrt(fmax(0., kallen_lambda(m0 * m0, m1 * m1, m2 * m2))) / m0;
}

// Two-body decay - isotropic
__device__ bool one_to_two_decay(vec4 p0, double m1, double m2, vec4& p1,
                                vec4& p2, double rho_1, double rho_2);

// Two-body decay - axis
__device__ bool one_to_two_decay(vec4 p0, double m1, double m2, vec4& p1,
                                vec4& p2, vec4 axis);

// N Body decay
__device__ bool one_to_n_decay(vec4 p0, const double* masses, int n,
                              vec4* momenta, event& ev);

#endif  // decay_kinematics_cuh_
