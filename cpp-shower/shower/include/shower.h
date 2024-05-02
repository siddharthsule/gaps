#ifndef shower_h_
#define shower_h_

#include "event.h"
#include "qcd.h"

// splitting function codes - only ff for now (removed zeroes)
// ------------------------------------------
const int sf_codes[] = {1,  2,  3,   4,   5,   11,  12,  13,
                        14, 15, 200, 301, 302, 303, 304, 305};

// shower class - not in gpu as kernels cannot be member functions
class shower {
 private:
  alpha_s as = alpha_s(mz, asmz);

 public:
  shower();  // t_c and asmax are preset in base

  // in cuda, we cannot point by reference, so we use
  // pointers (*) instead of references (&).

  // splitting functions
  double sf_value(double z, double y, int sf);
  double sf_estimate(double z, int sf);
  double sf_integral(double zm, double zp, int sf);
  double sf_generate_z(double zm, double zp, double rand, int sf);

  // utilities
  bool validate_splitting(int ij, const int sf);
  void sf_to_flavs(int sf, int* flavs);

  // select the winner emission
  void select_winner(event& ev, std::mt19937& gen);

  // kinematics
  void make_kinematics(vec4* kinematics, const double z, const double y,
                       const double phi, const vec4 pijt, const vec4 pkt);

  // colours
  void make_colours(event& ev, int* coli, int* colj, const int flavs[3],
                    const int colij[2], const int colk[2], const double r);

  // veto algorithm + perform the splitting
  void generate_splitting(event& ev, std::mt19937& gen);

  // run the shower
  void run(event& ev);
};

#endif  // shower_h_