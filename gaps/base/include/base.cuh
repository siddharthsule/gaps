#ifndef BASE_CUH_
#define BASE_CUH_

/**
 * The Base Class
 * --------------
 *
 * This file contains the neccessary includes and definitions that are used
 * throughout the program. This includes CUDA Libraries, Thrust Libraries,
 * C++ Libraries, and some global variables. Make changes here if you want to
 * change the global settings of the program!
 *
 * (but also be careful with what you change, as it may break the program...)
 */

// -----------------------------------------------------------------------------
// Import Libraries

// CUDA Libraries
#include <cuda_runtime.h>   // CUDA Runtime
#include <curand_kernel.h>  // CURAND Library

// Thrust Libraries
#include <thrust/device_vector.h>  // ALL EVENTS ON DEVICE

// C++ Libraries (Genreally Used)
#include <cmath>     // Math Functions
#include <cstdlib>   // SYS EXIT Command
#include <fstream>   // File I/O
#include <iostream>  // Standard I/O

// -----------------------------------------------------------------------------
// Program Settings - CAREFUL WITH CHANGES

// Debugging - only debug levels 0 and 1 (true or false)
const bool debug = false;

// Max Number of Partons, set to save memory
// at 10^6 Events:
//  50 works for all, but observables calc is slow
//  100 works for ME + PS, but not for Observables
//  30 is more than enough to do e+e- at 91.2 GeV
const int maxPartons = 30;

// LEP 91.2 settings
const double mz = 91.1876;
const double asmz = 0.118;

// Cutoff and its value of alpha_s (pre-calculated)
const double tC = 1.;
const double asmax = 0.440886;

// Number of Histogram Bins: Common for all Plots (for now...)
const int nBins = 100;
const int nBins2D = 10;  // 10x10 Grid

// Maximum Number of Events, beyond which program will be done in batches
const int maxEvents = 1000000;

// -----------------------------------------------------------------------------
// Common Functions

// Sync Device and Check for CUDA Errors
void syncGPUAndCheck(const char *operation);

// Debugging Function - Available in Kernels too!
__host__ __device__ void DEBUG_MSG(const char *msg);

// -----------------------------------------------------------------------------
#endif  // BASE_CUH_