#ifndef MATRIX_H_
#define MATRIX_H_

// Parton includes Base, which has the CUDA libraries
#include "parton.h"

// Matrix class
class Matrix {
 public:
  // Constructor
  Matrix(double alphas, double ecms = 91.2);

  // Leading Order Matrix Element Generation
  double ME2(int fl, double s, double t);

  // Generate a leading order point
  void GenerateLOPoint(Event &ev);

 private:
  double alphas, ecms, MZ2, GZ2, alpha, sin2tw, amin, ye, ze, ws;
};

#endif  // MATRIX_H_