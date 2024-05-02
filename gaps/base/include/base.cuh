#ifndef base_cuh_
#define base_cuh_

/**
 * the base class
 * --------------
 *
 * this file contains the neccessary includes and definitions that are used
 * throughout the program. this includes cuda libraries, thrust libraries,
 * c++ libraries, and some global variables. make changes here if you want to
 * change the global settings of the program!
 *
 * (but also be careful with what you change, as it may break the program...)
 */

// -----------------------------------------------------------------------------
// import libraries

// cuda libraries
#include <cuda_runtime.h>   // cuda runtime
#include <curand_kernel.h>  // curand library

// thrust libraries
#include <thrust/device_vector.h>  // all events on device

// c++ libraries (genreally used)
#include <cmath>     // math functions
#include <cstdlib>   // sys exit command
#include <fstream>   // file i/o
#include <iostream>  // standard i/o

// -----------------------------------------------------------------------------
// program settings - careful with changes

// debugging - only debug levels 0 and 1 (true or false)
const bool debug = false;

// max number of partons, set to save memory
// at 10^6 events:
//  50 works for all, but observables calc is slow
//  100 works for me + ps, but not for observables
//  30 is more than enough to do e+e- at 91.2 ge_v
const int max_partons = 30;

// lep 91.2 settings
const double mz = 91.1876;
const double asmz = 0.118;

// cutoff and its value of alpha_s (pre-calculated)
const double t_c = 1.;
const double asmax = 0.440886;

// number of histogram bins: common for all plots (for now...)
const int n_bins = 100;
const int n_bins2d = 100;  // nx_n grid

// maximum number of events, beyond which program will be done in batches
const int max_events = 1000000;

// -----------------------------------------------------------------------------
// common functions

// sync device and check for cuda errors
void sync_gpu_and_check(const char *operation);

// debugging function - available in kernels too!
__host__ __device__ void debug_msg(const char *msg);

#endif  // base_cuh_