#ifndef shower_cuh_
#define shower_cuh_

// qcd includes all the necessary headers
#include "qcd.cuh"

/**
 * a dipole shower on gpu
 * ----------------------
 *
 * note: kernel = cuda function, splitting function = qcd function
 *
 * this is the main result of the published work. it is a full implementation of
 * a dipole shower on the gpu. it is designed to be as fast as possible*, and
 * uses a number of tricks to achieve this. the main trick is to use a single
 * kernel to perform the entire shower, and to use a number of optimisations to
 * make the code as fast as possible.
 *
 * with the event object storing all the neccessary information and with the
 * fact that kernel's can't be member functions, the shower class has been
 * removed
 *
 * *: as possible as a second year ph_d student can make it ;)
 */

// prepare events for the shower
__global__ void prep_shower(event *events, int n);

// selecting the winner emission
__global__ void select_winner_split_func(event *events, int n);

// check the cutoff
__global__ void check_cutoff(event *events, int *d_completed, double cutoff,
                             int n);

// veto algorithm
__global__ void veto_alg(event *events, double *asval, bool *accept_emission,
                         int n);

// perform the splitting
__global__ void do_splitting(event *events, bool *accept_emission, int n);

// kinematics
__device__ void make_kinematics(vec4 *kinematics, const double z,
                                const double y, const double phi,
                                const vec4 pijt, const vec4 pkt);

// all tasks wrapped into a function
void run_shower(thrust::device_vector<event> &d_events);

#endif  // shower_cuh_
