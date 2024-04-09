#ifndef SHOWER_CUH_
#define SHOWER_CUH_

// qcd includes all the necessary headers
#include "qcd.cuh"

// Shower Kinematics and Colour in separate files
#include "shower_kinematics.cuh"
#include "colours.cuh"

/**
 * A Dipole Shower on GPU
 * ----------------------
 *
 * NOTE: Kernel = CUDA function, Splitting Function = QCD function
 *
 * This is the main result of the published work. It is a full implementation of
 * a dipole shower on the GPU. It is designed to be as fast as possible*, and
 * uses a number of tricks to achieve this. The main trick is to use a single
 * kernel to perform the entire shower, and to use a number of optimisations to
 * make the code as fast as possible.
 *
 * With the Event Object storing all the neccessary information and with the
 * fact that kernel's can't be member functions, the Shower Class has been
 * removed
 *
 * *: as possible as a second year PhD student can make it ;)
 */

// Initialise the curandStates
__global__ void initCurandStates(curandState *states, int N);

// Prepare Events for the Shower
__global__ void prepShower(Event *events, int N);

// Selecting the Winner Emission
__global__ void selectWinnerSplitFunc(Event *events, curandState *states,
                                      int N);

// Check the Cutoff
__global__ void checkCutoff(Event *events, int *d_completed, double cutoff,
                            int N);

// Veto Algorithm
__global__ void vetoAlg(Event *events, curandState *states, int N);

// Perform the Splitting
__global__ void doSplitting(Event *events, curandState *states, int N);

// Kinematics
__device__ void MakeKinematics(Vec4 *kinematics, const double z, const double y,
                               const double phi, const Vec4 pijt,
                               const Vec4 pkt);

// All tasks wrapped into a function
void runShower(thrust::device_vector<Event> &d_events);

#endif  // SHOWER_CUH_
