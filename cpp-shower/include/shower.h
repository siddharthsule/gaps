#ifndef SHOWER_H_
#define SHOWER_H_

#include "parton.h"
#include "qcd.h"

/**
 * Removal of Splitting Function Class
 * -----------------------------------
 *
 * As we are doing flavourless quarks, to make the GPU implementation easier,
 * we will remove the Splitting Function Class. This class only contains const
 * functions, so we can just use the functions directly in the Shower class.
 *
 * For more detailed implementations, the Splitting Function Class can be
 * brought back, but we really dont need it for now.
 */

// Shower Class - NOT IN GPU as kernels cannot be member functions
class Shower {
 public:
  Shower(); // tC and asmax are preset in base

  // In CUDA, we cannot point by reference, so we use
  // pointers (*) instead of references (&).
  void SelectWinner(Event& ev, std::mt19937& gen);
  void MakeKinematics(Vec4* kinematics, const double z, const double y,
                      const double phi, const Vec4 pijt, const Vec4 pkt);
  void MakeColours(Event& ev, int* coli, int* colj, const int flavs[3],
                   const int colij[2], const int colk[2], const int r);
  void GenerateSplitting(Event& ev, std::mt19937& gen);
  void Run(Event& ev, double t);

 private:
  AlphaS as = AlphaS(mz, asmz);
};

#endif  // SHOWER_H_