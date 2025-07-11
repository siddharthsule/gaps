#ifndef shower_cuh_
#define shower_cuh_

#include "qcd.cuh"

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
 public:
  double t_c;
  double as_max;
  double j_max;

 public:
  // constructor
  __device__ void setup(double root_s, double t_c, double as_max);

  // splitting functions
  __device__ double sf_value(double z, double y, int sf) const;
  __device__ double sf_estimate(double z, int sf) const;
  __device__ double sf_integral(double zm, double zp, int sf) const;
  __device__ double sf_generate_z(double zm, double zp, double rand,
                                  int sf) const;

  // utility functions to switch between codes and splitting functions
  __device__ int get_splitting_case(int sf) const;
  __device__ int get_splitting_type(int sf) const;
  __device__ int is_emitter_antiparticle(int sf) const;
  __device__ int get_splitting_flavour(int sf) const;

  // splitting functions for the shower
  __device__ void get_possible_splittings(int ij, int *splittings) const;
  __device__ bool validate_splitting(int ij, int sf, bool emt_init,
                                     bool spc_init) const;
  __device__ void sf_to_flavs(int sf, int *flavs) const;
  __device__ void sf_to_text(int sf, char *text) const;

  // kinematics
  __device__ void make_kinematics(vec4 *kinematics, const double z,
                                  const double y, const double phi,
                                  const vec4 pijt, const vec4 pkt,
                                  int sf) const;

  // colours
  __device__ void make_colours(int current_col, int sf, int flavs[3],
                               int colij[2], int colk[2], int *coli, int *colj,
                               double r) const;
};

// -----------------------------------------------------------------------------
// Kernels and Wrappers

__global__ void shower_setup_kernel(shower *sh, double as_max, double root_s);

__global__ void prep_shower(event *events, bool nlo_matching, int n);

__global__ void select_winner_split_func(shower *shower, event *events, int n,
                                         double *winner);

__global__ void check_cutoff(event *events, int *d_completed, double cutoff,
                             int n);

__global__ void veto_alg(shower *shower, alpha_s *as, event *events, int n,
                         double *xf_a, double *xf_b, bool *accept_emission,
                         double *winner, int *d_evaluations,
                         int *d_overestimate_error);

__global__ void do_splitting(shower *shower, event *events, int n,
                             bool *accept_emission, double *winner);

__global__ void check_too_many_particles(event *events,
                                         int *d_too_many_particles, int n);

// all tasks wrapped into a function
void run_shower(thrust::device_vector<event> &dv_events, double root_s,
                bool nlo_matching, bool do_partition, double t_c, double asmz,
                int n_emissions_max, int blocks, int threads);

#endif  // shower_cuh_