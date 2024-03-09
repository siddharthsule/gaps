#ifndef SHOWER_H_
#define SHOWER_H_

#include "parton.h"
#include "qcd.h"

// Shower Class - NOT IN GPU as kernels cannot be member functions
class Shower {
 private:
  AlphaS as = AlphaS(mz, asmz);

 public:
  Shower();  // tC and asmax are preset in base

  // In CUDA, we cannot point by reference, so we use
  // pointers (*) instead of references (&).

  // Select the Winner Emission
  void SelectWinner(Event& ev, std::mt19937& gen);

  // Kinematics
  void MakeKinematics(Vec4* kinematics, const double z, const double y,
                      const double phi, const Vec4 pijt, const Vec4 pkt);

  // Colours
  void MakeColours(Event& ev, int* coli, int* colj, const int flavs[3],
                   const int colij[2], const int colk[2], const int r);

  // Veto Algorithm + Perform the Splitting
  void GenerateSplitting(Event& ev, std::mt19937& gen);

  // Run the Shower
  void Run(Event& ev);
};

#endif  // SHOWER_H_