#ifndef base_h_
#define base_h_

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

#include <cmath>     // math functions
#include <fstream>   // file i/o
#include <iostream>  // standard i/o
// #include <random>    // random number generation
#include <vector>  // vector

// -----------------------------------------------------------------------------
// program settings - careful with changes

// RNG Settings - Linear Congruential Generator
const unsigned long lcg_a = 1664525;
const unsigned long lcg_c = 1013904223;
const unsigned long lcg_m = 4294967296;

// max number of partons, set to save memory
// at 10^6 events:
//  50 works for all, but observables calc is slow
//  100 works for me + ps, but not for observables
//  30 is more than enough to do e+e- at 91.2 GeV
const int max_partons = 30;

// lep 91.2 settings
const double mz = 91.1876;
const double asmz = 0.118;

// cutoff and its value of alpha_s (pre-calculated)
const double t_c = 1.;
const double asmax = 0.440886;

// number of histogram bins: common for all plots (for now...)
const int n_bins = 100;
const int n_bins2d = 100;  // 10x10 grid

// maximum number of events, beyond which program will be done in batches
const int max_events = 1000000;

#endif  // base_h_