#ifndef base_cuh_
#define base_cuh_

/**
 * the base file
 * -------------
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

// c++ libraries (generally used)
#include <cmath>     // math functions
#include <cstdlib>   // standard library
#include <fstream>   // file i/o
#include <iostream>  // standard i/o

// cuda libraries
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/partition.h>

// -----------------------------------------------------------------------------
// program settings - careful with changes

// maximum number of events the GPU can handle
const int max_events = 1000000;

// max number of particles, set to 50 to save memory
const int max_particles = 50;

// RNG Settings - Linear Congruential Generator
const unsigned long lcg_a = 1664525;
const unsigned long lcg_c = 1013904223;
const unsigned long lcg_m = 4294967296;

// QCD: number of colors and flavors
const double k_nc = 3.;
const int k_nf = 5;

// Z Boson: Mass, Width
const double mz = 91.1876;
const double gz = 2.4952;

// XS: GeV^-2 to pb conversion
const double GeV_minus_2_to_pb = 3.89379656e8;

// Observables: number of histogram bins
const int n_bins = 100;

// -----------------------------------------------------------------------------
// common functions and debugging

// sync device and check for cuda errors
void sync_gpu_and_check(const char *operation);

// debugging - only debug levels 0 and 1 (true or false)
const bool debug = false;

// debugging function - available in kernels too!
__host__ __device__ void debug_msg(const char *msg);

#endif  // base_cuh_