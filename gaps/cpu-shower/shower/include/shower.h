#ifndef shower_h_
#define shower_h_

// event and qcd includes all the necessary headers
#include "event.h"
#include "qcd.h"

class shower {
  /**
   * @class shower
   * @brief the dipole shower
   *
   * this is the main result of the published work. it is a full implementation
   * of a dipole shower on the gpu. it is designed to be as fast as possible*,
   * and uses a number of tricks to achieve this. the main trick is to use a
   * single kernel to perform the entire shower, and to use a number of
   * optimisations to make the code as fast as possible.
   *
   * with the event object storing all the neccessary information and with the
   * fact that kernel's can't be member functions, the shower class has been
   * removed
   *
   * *: as possible as a second year phd student can make it ;)
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables
  alpha_s as;

 public:
  double t_c;
  double as_max;
  double j_max;

  int n_emissions_max;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  shower(double root_s, double t_c, double asmz, int n_emissions_max)
      : t_c(t_c), n_emissions_max(n_emissions_max), as(mz, asmz) {
    as_max = as(t_c);
    j_max = 1.;
  }

  // ---------------------------------------------------------------------------
  // member functions

  // splitting functions
  double sf_value(double z, double y, int sf) const;
  double sf_estimate(double z, int sf) const;
  double sf_integral(double zm, double zp, int sf) const;
  double sf_generate_z(double zm, double zp, double rand, int sf) const;

  // utility functions to switch between codes and splitting functions
  int get_splitting_case(int sf) const;
  int get_splitting_type(int sf) const;
  int is_emitter_antiparticle(int sf) const;
  int get_splitting_flavour(int sf) const;

  // splitting functions for the shower
  void get_possible_splittings(int ij, int* splittings) const;
  bool validate_splitting(int ij, const int sf, bool emt_init,
                          bool spc_init) const;
  void sf_to_flavs(int sf, int* flavs) const;
  void sf_to_text(int sf, char* text) const;

  // select the winner emission
  void select_winner(event& ev, double* winner) const;

  // kinematics
  void make_kinematics(vec4* kinematics, const double z, const double y,
                       const double phi, const vec4 pijt, const vec4 pkt,
                       int sf) const;

  // colours
  void make_colours(int current_col, int sf, int flavs[3], int colij[2],
                    int colk[2], int* coli, int* colj, double r) const;

  // veto algorithm + perform the splitting
  void generate_splitting(event& ev);

  // run the shower
  void run(event& ev, bool nlo_matching);
};

#endif  // shower_h_