#ifndef event_h_
#define event_h_

#include "parton.h"

// event class
// built to contain the partons and the dxs as one accesible object
// in future, can use to store thrust, log10y23 to parallelise those
class event {
 private:
  // temporary solution - allows a limited number of partons
  // better solution would be to use a dynamic array, but not gpu friendly
  parton partons[max_partons];

  // random seed ---------------------------------------------------------------

  unsigned long seed = 0;
  double rand = 0.;

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
  event() {}

  // getters -------------------------------------------------------------------

  // get seed and random number
  unsigned long get_seed() const { return seed; }
  double get_rand() const { return rand; }

  // access partons in the event
  parton get_parton(int i) const { return partons[i]; }
  int get_size() const { return n_hard + n_emission; }
  int get_hard() const { return n_hard; }
  int get_emissions() const { return n_emission; }
  int get_parton_size() const {
    return (n_hard + n_emission) - 2;
  }  // -2: e+, e-

  // get differential cross section
  double get_dxs() const { return dxs; }

  // get shower params
  double get_shower_t() const { return shower_t; }
  double get_shower_y() const { return shower_y; }
  double get_shower_z() const { return shower_z; }
  int get_shower_c() const { return shower_c; }

  // get winner emission
  int get_win_sf() const { return win_sf; }
  int get_win_dipole(int i) const { return win_dipole[i]; }
  double get_win_param(int i) const { return win_params[i]; }

  // get analysis variables
  bool get_validity() const { return validity; }

  double get_y23() const { return y23; }
  double get_y34() const { return y34; }
  double get_y45() const { return y45; }
  double get_y56() const { return y56; }
  double get_thr() const { return thr; }
  double get_hjm() const { return hjm; }
  double get_ljm() const { return ljm; }
  double get_wjb() const { return wjb; }
  double get_njb() const { return njb; }

  vec4 get_t_axis() const { return t_axis; }

  double get_dalitz(int i) const { return dalitz[i]; }

  // setters -------------------------------------------------------------------

  // set random seed and random number
  void set_seed(unsigned long seed) { this->seed = seed; }
  void set_rand(double rand) { this->rand = rand; }

  // add / replace parton
  void set_parton(int i, parton parton) { partons[i] = parton; }

  // not used in me [host]
  void set_parton_pid(int i, int pid) { partons[i].set_pid(pid); }
  void set_parton_mom(int i, vec4 mom) { partons[i].set_mom(mom); }
  void set_parton_col(int i, int col) { partons[i].set_col(col); }
  void set_parton_anti_col(int i, int anticol) {
    partons[i].set_anti_col(anticol);
  }

  // set differential cross section and n_hard
  void set_dxs(double dxs) { this->dxs = dxs; }
  void set_hard(int n_hard) { this->n_hard = n_hard; }

  // adjust and increment number of emissions
  void set_emissions(int n_emission) { this->n_emission = n_emission; }
  void increment_emissions() { n_emission++; }

  // set shower params
  void set_shower_t(double shower_t) { this->shower_t = shower_t; }
  void set_shower_y(double shower_y) { this->shower_y = shower_y; }
  void set_shower_z(double shower_z) { this->shower_z = shower_z; }

  void set_shower_c(int shower_c) { this->shower_c = shower_c; }
  void increment_shower_c() { shower_c++; }

  // set winner emission
  void set_win_sf(int win_sf) { this->win_sf = win_sf; }
  void set_win_dipole(int i, int win_dipole) {
    this->win_dipole[i] = win_dipole;
  }
  void set_win_param(int i, double win_params) {
    this->win_params[i] = win_params;
  }

  // set analysis variables
  void set_validity(bool validity) { this->validity = validity; }

  void set_y23(double y23) { this->y23 = y23; }
  void set_y34(double y34) { this->y34 = y34; }
  void set_y45(double y45) { this->y45 = y45; }
  void set_y56(double y56) { this->y56 = y56; }
  void set_thr(double thr) { this->thr = thr; }
  void set_hjm(double hjm) { this->hjm = hjm; }
  void set_ljm(double ljm) { this->ljm = ljm; }
  void set_wjb(double wjb) { this->wjb = wjb; }
  void set_njb(double njb) { this->njb = njb; }

  void set_t_axis(vec4 t_axis) { this->t_axis = t_axis; }

  void set_dalitz(double x1, double x2) {
    dalitz[0] = x1;
    dalitz[1] = x2;
  }

  // member functions ----------------------------------------------------------

  // validation of result data
  bool validate() {
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

    bool pcheck = (std::abs(psum[0]) < 1e-12 && std::abs(psum[1]) < 1e-12 &&
                   std::abs(psum[2]) < 1e-12 && std::abs(psum[3]) < 1e-12);
    if (!pcheck) {
      std::cout << psum << std::endl;
    }

    bool ccheck = true;
    for (int i = 0; i < max_partons - 1; i++) {
      if (csum[i] != 0) {
        std::cout << "colour " << i << " is not conserved." << std::endl;
        ccheck = false;
        break;
      }
    }

    return pcheck && ccheck;
  }

  void print_info() const {
    std::cout << "event information:\n";
    std::cout << "dxs: " << get_dxs() << "\n";
    std::cout << "number of emissions: " << get_emissions() << "\n";
    std::cout << "shower t: " << get_shower_t() << "\n";
    std::cout << "shower y: " << get_shower_y() << "\n";
    std::cout << "shower z: " << get_shower_z() << "\n";
    std::cout << "shower c: " << get_shower_c() << "\n";
    std::cout << "winner sf: " << get_win_sf() << "\n";
    std::cout << "winner dipole 1: " << get_win_dipole(0) << "\n";
    std::cout << "winner dipole 2: " << get_win_dipole(1) << "\n";
    std::cout << "winner params 1: " << get_win_param(0) << "\n";
    std::cout << "winner params 2: " << get_win_param(1) << "\n";
    std::cout << "partons:\n";
    for (int i = 0; i < get_size(); i++) {
      parton parton = get_parton(i);
      std::cout << "  parton " << i << ":\n";
      std::cout << "    pid: " << parton.get_pid() << "\n";
      std::cout << "    mom: " << parton.get_mom() << "\n";
      std::cout << "    col: " << parton.get_col() << "\n";
      std::cout << "    anti_col: " << parton.get_anti_col() << "\n";
    }
  }
};

#endif  // event_h_