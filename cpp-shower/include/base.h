#ifndef BASE_H_
#define BASE_H_

#include <cmath> // Math Functions
#include <vector> // Vector
#include <random> // Random Number Generation
#include <fstream> // File I/O
#include <iostream> // Standard I/O

// Max Number of Partons, set to save memory
const int maxPartons = 30;

// LEP 91.2 settings
const double mz = 91.1876;
const double asmz = 0.118;

// Cutoff and its value of alpha_s (pre-calculated)
const double tC = 1.0;
const double asmax = 0.440886;

// Number of Histogram Bins: Common for all Plots (for now...)
const int nBins = 100;

#endif  // BASE_H_