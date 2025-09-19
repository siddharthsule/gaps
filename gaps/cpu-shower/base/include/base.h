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

// bring math functions into global namespace for convenience
using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::ceil;
using std::cos;
using std::exp;
using std::fabs;
using std::floor;
using std::isinf;
using std::isnan;
using std::log;
using std::max;
using std::min;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::tan;

// LHAPDF
#include <LHAPDF/LHAPDF.h>

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

// Charm and Bottom quark masses
const double mc = 1.30;
const double mb = 4.75;

// Z Boson: Mass, Width
const double mz = 91.1876;
const double gz = 2.4952;

// PDF: Limits (CT14lo)
const double pdf_x_min = 1e-9;
const double pdf_x_max = 1.;
const double pdf_q_min = 1.295;
const double pdf_q_max = 1e5;

// ME: Z mass cuts for qq -> ee
const double mz_cut_a = 60.;
const double mz_cut_b = 120.;

// XS: GeV^-2 to pb conversion
const double GeV_minus_2_to_pb = 3.89379656e8;

// Observables: max number of histogram bins
const int max_bins = 100;

// Gen Kt: power and R values
const double power = -1.;  // 1. = kt, 0. = cambridge, -1. = anti-kt
const double R = 0.4;

#endif  // base_h_