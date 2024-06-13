#ifndef shower_h_
#define shower_h_

#include "event.h"
#include "qcd.h"

/**
 * a dipole shower on gpu
 * ----------------------
 *
 * note: kernel = cuda function, splitting function = qcd function
 *
 * this is the main result of the published work. it is a full implementation of
 * a dipole shower on the gpu. it is designed to be as fast as possible*, and
 * uses a number of tricks to achieve this. the main trick is to use a single
 * kernel to perform the entire shower, and to use a number of optimisations to
 * make the code as fast as possible.
 *
 * with the event object storing all the neccessary information and with the
 * fact that kernel's can't be member functions, the shower class has been
 * removed
 *
 * *: as possible as a second year ph_d student can make it ;)
 */

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
  void select_winner(event& ev);

  // kinematics
  void make_kinematics(vec4* kinematics, const double z, const double y,
                       const double phi, const vec4 pijt, const vec4 pkt);

  // colours
  void make_colours(event& ev, int* coli, int* colj, const int flavs[3],
                    const int colij[2], const int colk[2], const double r);

  // veto algorithm + perform the splitting
  void generate_splitting(event& ev);

  // run the shower
  void run(event& ev);
};

#endif  // shower_h_