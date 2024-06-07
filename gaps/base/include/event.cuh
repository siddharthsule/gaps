#ifndef event_cuh_
#define event_cuh_

#include "parton.cuh"
#include "prng.cuh"

/**
 * the event class
 * ---------------
 *
 * this is the backbone of the program. an event contains the differential cross
 * section calculated from the me, the hard partons and the showered partons.
 * it also contains the shower parameters which are constantly updated during
 * the showering process. the event also contains the analysis variables, which
 * are calculated after the showering process is complete.
 */

class event {
 private:
  // random seed ---------------------------------------------------------------

  unsigned long seed = 0;

  // partons -------------------------------------------------------------------

  // temporary solution - allows a limited number of partons
  // better solution would be to use a dynamic array, but not gpu friendly
  parton partons[max_partons];

  // me params -----------------------------------------------------------------

  double dxs = 0.;  // differential cross section
  int n_hard = 0;   // number of hard partons
  // int n_initial = 0;    // number of initial partons (prep for isr)
  // int n_non_parton = 0;  // number of non-parton partons (prep for isr)

  // shower params -------------------------------------------------------------

  int n_emission = 0;    // number of emissions
  double shower_t = 0.;  // evolution and splitting variables
  double shower_z = 0.;
  double shower_y = 0.;
  int shower_c = 0;  // colour counter

  // selecting winner emission - defaults values which represent no winner
  int win_sf = 16;
  int win_dipole[2] = {-1, -1};
  double win_params[2] = {0., 0.};

  // end shower flag
  bool end_shower = false;  // shower end flag - used if t < 1 ge_v

  // analysis variables --------------------------------------------------------

  // event validity - momentum and colour conservation
  bool validity = true;

  // jet rates using the durham algorithm
  double y23 = -50., y34 = -50., y45 = -50., y56 = -50.;

  // event shape variables - thrust, jet masses and broadenings
  double thr = -50., hjm = -50., ljm = -50., wjb = -50., njb = -50.;
  vec4 t_axis = vec4();

  // dalitz plot
  double dalitz[2] = {-50., -50.};

 public:
  // constructor ---------------------------------------------------------------

  // empty, so that we can build our me, ps onto it
  __device__ event() {}

  // getters -------------------------------------------------------------------

  // random seed
  __device__ unsigned long get_seed() const { return seed; }

  // access partons in the event
  __device__ parton get_parton(int i) const { return partons[i]; }
  __device__ int get_size() const { return n_hard + n_emission; }
  __device__ int get_hard() const { return n_hard; }
  __device__ int get_emissions() const { return n_emission; }
  __device__ int get_parton_size() const {
    // -2: e+, e-
    return (n_hard + n_emission) - 2;
  }

  // get differential cross section
  __device__ double get_dxs() const { return dxs; }

  // get shower params
  __device__ double get_shower_t() const { return shower_t; }
  __device__ double get_shower_z() const { return shower_z; }
  __device__ double get_shower_y() const { return shower_y; }
  __device__ int get_shower_c() const { return shower_c; }

  __device__ int get_win_sf() const { return win_sf; }
  __device__ int get_win_dipole(int i) const { return win_dipole[i]; }
  __device__ double get_win_param(int i) const { return win_params[i]; }

  __device__ bool get_end_shower() const { return end_shower; }

  // analysis getters
  __device__ bool get_validity() const { return validity; }

  __device__ double get_y23() const { return y23; }
  __device__ double get_y34() const { return y34; }
  __device__ double get_y45() const { return y45; }
  __device__ double get_y56() const { return y56; }
  __device__ double get_thr() const { return thr; }
  __device__ double get_hjm() const { return hjm; }
  __device__ double get_ljm() const { return ljm; }
  __device__ double get_wjb() const { return wjb; }
  __device__ double get_njb() const { return njb; }

  __device__ vec4 get_t_axis() const { return t_axis; }

  __device__ double get_dalitz(int i) const { return dalitz[i]; }

  // setters -------------------------------------------------------------------

  // set random seed
  __device__ void set_seed(unsigned long seed) { this->seed = seed; }

  // add / replace parton
  __device__ void set_parton(int i, parton parton) { partons[i] = parton; }

  // set parton data
  __device__ void set_parton_pid(int i, int pid) { partons[i].set_pid(pid); }
  __device__ void set_parton_mom(int i, vec4 mom) { partons[i].set_mom(mom); }
  __device__ void set_parton_col(int i, int col) { partons[i].set_col(col); }
  __device__ void set_parton_anti_col(int i, int anticol) {
    partons[i].set_anti_col(anticol);
  }

  // set differential cross section and n_hard
  __device__ void set_dxs(double dxs) { this->dxs = dxs; }
  __device__ void set_hard(int n_hard) { this->n_hard = n_hard; }

  // adjust and increment number of emissions
  __device__ void set_emissions(int n_emission) {
    this->n_emission = n_emission;
  }
  __device__ void increment_emissions() { n_emission++; }

  // set shower params
  __device__ void set_shower_t(double shower_t) { this->shower_t = shower_t; }
  __device__ void set_shower_z(double shower_z) { this->shower_z = shower_z; }
  __device__ void set_shower_y(double shower_y) { this->shower_y = shower_y; }

  __device__ void set_shower_c(int shower_c) { this->shower_c = shower_c; }
  __device__ void increment_shower_c() { shower_c++; }

  __device__ void set_win_sf(int win_sf) { this->win_sf = win_sf; }
  __device__ void set_win_dipole(int i, int win_parton) {
    this->win_dipole[i] = win_parton;
  }
  __device__ void set_win_param(int i, double win_param) {
    this->win_params[i] = win_param;
  }

  __device__ void set_end_shower(bool end_shower) {
    this->end_shower = end_shower;
  }

  // set analysis variables
  __device__ void set_validity(bool validity) { this->validity = validity; }

  __device__ void set_y23(double y23) { this->y23 = y23; }
  __device__ void set_y34(double y34) { this->y34 = y34; }
  __device__ void set_y45(double y45) { this->y45 = y45; }
  __device__ void set_y56(double y56) { this->y56 = y56; }
  __device__ void set_thr(double thr) { this->thr = thr; }
  __device__ void set_hjm(double hjm) { this->hjm = hjm; }
  __device__ void set_ljm(double ljm) { this->ljm = ljm; }
  __device__ void set_wjb(double wjb) { this->wjb = wjb; }
  __device__ void set_njb(double njb) { this->njb = njb; }

  __device__ void set_t_axis(vec4 t_axis) { this->t_axis = t_axis; }

  __device__ void set_dalitz(double x1, double x2) {
    dalitz[0] = x1;
    dalitz[1] = x2;
  }

  // member functions ----------------------------------------------------------

  // validate the event - check momentum and colour conservation
  __device__ bool validate() {
    vec4 psum = vec4();

    // n colours = n partons - 1
    int csum[max_partons - 1] = {0};

    for (int i = 0; i < get_size(); i++) {
      parton p = get_parton(i);

      vec4 pmom = p.get_mom();
      int pcol = p.get_col();
      int p_anti_col = p.get_anti_col();

      psum = psum + pmom;

      if (pcol > 0) {
        csum[pcol] += 1;
      }

      if (p_anti_col > 0) {
        csum[p_anti_col] -= 1;
      }
    }

    bool pcheck = (abs(psum[0]) < 1e-12 && abs(psum[1]) < 1e-12 &&
                   abs(psum[2]) < 1e-12 && abs(psum[3]) < 1e-12);

    /* // no need to print for gpu, it counts number of invalid events
    if (!pcheck) {
      printf("%f %f %f %f\n", psum[0], psum[1], psum[2], psum[3]);
    }
    */

    bool ccheck = true;
    for (int i = 0; i < max_partons - 1; i++) {
      if (csum[i] != 0) {
        // printf("colour %d is not conserved.\n", i);
        ccheck = false;
        break;
      }
    }

    return pcheck && ccheck;
  }

  __device__ void print_info() const {
    printf("event information:\n");
    printf("dxs: %f\n", get_dxs());
    printf("number of emissions: %d\n", get_emissions());
    printf("shower parameters:\n");
    printf("  t: %f\n", get_shower_t());
    printf("  y: %f\n", get_shower_y());
    printf("  z: %f\n", get_shower_z());
    printf("  c: %d\n", get_shower_c());
    printf("shower winner:\n");
    printf("  kernel number: %d\n", get_win_sf());
    printf("  partons: [%d, %d]\n", get_win_dipole(0), get_win_dipole(1));
    printf("  params: [%f, %f]\n", get_win_param(0), get_win_param(1));
    printf("partons:\n");
    for (int i = 0; i < get_size(); i++) {
      parton parton = get_parton(i);
      printf("  parton %d:\n", i);
      printf("    pid: %d\n", parton.get_pid());
      printf("    mom: %f\n", parton.get_mom().p());
      printf("    col: %d\n", parton.get_col());
      printf("    anti_col: %d\n", parton.get_anti_col());
    }
  }

  // Generate a random seed and random number based on the Event's seed
  __device__ double gen_random() { return generate_lcg(seed); }
};

#endif  // event_cuh_