#ifndef base_h_
#define base_h_

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
#include <vector>    // vector container

// -----------------------------------------------------------------------------
// program settings - careful with changes

// maximum number of events the GPU can handle
const int max_events = 1000000;

// max number of particles, set to 50 to save memory
const int max_particles = 50;

// RNG: Linear Congruential Generator
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

#endif  // base_h_