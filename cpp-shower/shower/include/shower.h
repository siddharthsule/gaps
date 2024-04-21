#ifndef SHOWER_H_
#define SHOWER_H_

#include "event.h"
#include "qcd.h"

// Splitting Function Codes - Only FF for now (Removed Zeroes)
// ------------------------------------------
const int sfCodes[] = {1,  2,  3,   4,   5,   11,  12,  13,
                       14, 15, 200, 301, 302, 303, 304, 305};

// Shower Class - NOT IN GPU as kernels cannot be member functions
class Shower {
 private:
  AlphaS as = AlphaS(mz, asmz);

 public:
  Shower();  // tC and asmax are preset in base

  // In CUDA, we cannot point by reference, so we use
  // pointers (*) instead of references (&).

  // Splitting Functions
  double sfValue(double z, double y, int sf);
  double sfEstimate(double z, int sf);
  double sfIntegral(double zm, double zp, int sf);
  double sfGenerateZ(double zm, double zp, double rand, int sf);

  // Utilities
  bool validateSplitting(int ij, const int sf);
  void sfToFlavs(int sf, int* flavs);

  // Select the Winner Emission
  void SelectWinner(Event& ev, std::mt19937& gen);

  // Kinematics
  void MakeKinematics(Vec4* kinematics, const double z, const double y,
                      const double phi, const Vec4 pijt, const Vec4 pkt);

  // Colours
  void MakeColours(Event& ev, int* coli, int* colj, const int flavs[3],
                   const int colij[2], const int colk[2], const double r);

  // Veto Algorithm + Perform the Splitting
  void GenerateSplitting(Event& ev, std::mt19937& gen);

  // Run the Shower
  void Run(Event& ev, int seed = std::rand());
};

#endif  // SHOWER_H_