#ifndef event_h_
#define event_h_

#include "particle.h"
#include "prng.h"

class event {
  /**
   * @class event
   * @brief The event record.
   *
   * This is the backbone of the program. an event contains the differential
   * cross section calculated from the me, the hard particles and the showered
   * partons.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  // Event ID
  int id = 0;

  // random seed for PRNG, used to keep identical results to the GPU Shower
  unsigned long seed = 0;

  // Allows a limited number of particles (adjust max_particles in base.cuh)
  particle particles[max_particles];

  // Counters
  int n_hard = 0;      // Hard particles - Generated from the Matrix Element
  int n_emission = 0;  // Emissions - Generated from the Shower

  // event validity - momentum and colour conservation
  bool validity = true;

  // differential cross section, used as the event weight
  double me2 = 0.;
  double dxs = 0.;

  /**
   * I had thought of taking these out and putting them in the shower class, but
   * these params are used and updated for the whole shower and unique to each
   * event, it just made sense to me to keep them here. Each to their own!
   *
   * Alternative would be to have a struct for the shower params and pass it
   * around with the event. Completely fine alternative, but keeping them here
   * means that the event is self-contained and can be passed around without
   * needing to pass the shower params.
   */

  // shower parameters
  double shower_t = 0.;  // evolution variable
  int shower_c = 0;      // colour counter

  // shower end flag - used if t < t_cutoff
  bool end_shower = false;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  // empty, so that we can build our me, ps onto it
  event() {}

  // ---------------------------------------------------------------------------
  // getters

  int get_id() const {
    /**
     * @brief get the event id
     */

    return id;
  }

  unsigned long get_seed() const {
    /**
     * @brief get the event seed
     */

    return seed;
  }

  particle get_particle(int i) const {
    /**
     * @fn get_particle
     * @brief get the particle at index i of the event record
     *
     * @param i index of the particle
     */

    return particles[i];
  }

  int get_size() const {
    /**
     * @brief get the number of particles in the event record, equivalent to the
     * number of hard particles + the number of emissions
     */

    return n_hard + n_emission;
  }

  int get_hard() const {
    /**
     * @brief get the number of hard particles in the event record
     */

    return n_hard;
  }

  int get_emissions() const {
    /**
     * @brief get the number of emissions in the event record
     */

    return n_emission;
  }

  double get_me2() const {
    /**
     * @brief get the squared matrix element of the event
     */

    return me2;
  }

  double get_dxs() const {
    /**
     * @brief get the differential cross section of the event
     */

    return dxs;
  }

  double get_shower_t() const {
    /**
     * @brief get the shower evolution variable
     */

    return shower_t;
  }

  int get_shower_c() const {
    /**
     * @brief get the shower colour counter
     */

    return shower_c;
  }

  bool has_shower_ended() const {
    /**
     * @brief check if the shower has ended
     */

    return end_shower;
  }

  bool get_validity() const {
    /**
     * @brief get the validity of the event
     */

    return validity;
  }

  // ---------------------------------------------------------------------------
  // setters

  void set_id(int id) {
    /**
     * @brief set the event id
     *
     * @param id the event id
     */

    this->id = id;
  }

  void set_seed(unsigned long seed) {
    /**
     * @brief set the event seed
     *
     * @param seed the event seed
     */

    this->seed = seed;
  }

  void set_particle(int i, particle particle) {
    /**
     * @brief set the particle at index i of the event record. It completely
     * writes the (hopefully empty) particle at index i with the new particle.
     *
     * @param i index of the particle
     * @param particle the particle to set
     */

    particles[i] = particle;
  }

  /**
   * Setting particle Parameters through the Event?
   * --------------------------------------------
   *
   * I'm not sure if this is the best way to do this. I'm setting the particle
   * parameters through the event. This is because addresses on the GPU don't
   * work the same way as on the CPU. I can't pass the particle address to the
   * kernel and update the particle directly, As I am passing the event to the
   * kernel anyway, I can update the particle through the event
   *
   * Maybe in the future, someone can come up with a better way to do this.
   */

  void set_particle_pid(int i, int pid) {
    /**
     * @brief set the particle pid at index i of the event record
     *
     * @param i index of the particle
     * @param pid the particle pid
     */

    particles[i].set_pid(pid);
  }

  void set_particle_mom(int i, vec4 mom) {
    /**
     * @brief set the particle momentum at index i of the event record
     *
     * @param i index of the particle
     * @param mom the particle momentum
     */

    particles[i].set_mom(mom);
  }

  void set_particle_col(int i, int col) {
    /**
     * @brief set the particle colour at index i of the event record
     *
     * @param i index of the particle
     * @param col the particle colour
     */

    particles[i].set_col(col);
  }

  void set_particle_acol(int i, int acol) {
    /**
     * @brief set the particle anti-colour at index i of the event record
     *
     * @param i index of the particle
     * @param acol the particle anti-colour
     */

    particles[i].set_acol(acol);
  }

  void set_particle_eta(int i, double eta) {
    /**
     * @brief set the particle eta at index i of the event record
     *
     * @param i index of the particle
     * @param eta the particle eta
     */

    particles[i].set_eta(eta);
  }

  void set_hard(int n_hard) {
    /**
     * @brief set the number of hard particles in the event record
     *
     * @param n_hard the number of hard particles
     */

    this->n_hard = n_hard;
  }

  void set_emissions(int n_emission) {
    /**
     * @brief set the number of emissions in the event record
     *
     * @param n_emission the number of emissions
     */

    this->n_emission = n_emission;
  }

  void set_me2(double me2) {
    /**
     * @brief set the squared matrix element of the event
     *
     * @param me2 the squared matrix element
     */

    this->me2 = me2;
  }

  void set_dxs(double dxs) {
    /**
     * @brief set the differential cross section of the event
     *
     * @param dxs the differential cross section
     */

    this->dxs = dxs;
  }

  void set_shower_t(double shower_t) {
    /**
     * @brief set the shower evolution variable
     *
     * @param shower_t the shower evolution variable
     */

    this->shower_t = shower_t;
  }

  void set_shower_c(int shower_c) {
    /**
     * @brief set the shower colour counter
     *
     * @param shower_c the shower colour counter
     */

    this->shower_c = shower_c;
  }

  void increment_emissions() {
    /**
     * @brief increment the number of emissions and the shower colour counter
     */

    n_emission++;
    shower_c++;
  }

  void shower_has_ended(bool end_shower) {
    /**
     * @brief set the end shower flag
     *
     * @param end_shower the end shower flag
     */

    this->end_shower = end_shower;
  }

  void set_validity(bool validity) {
    /**
     * @brief set the validity of the event
     *
     * @param validity the validity of the event
     */

    this->validity = validity;
  }

  // ---------------------------------------------------------------------------
  // member functions

  void print_info() const {
    /**
     * @brief print the event information
     */

    printf("event number %d\n", get_id());
    printf("PRNG seed: %lu\n", get_seed());
    printf("\n");
    printf("dxs: %f\n", get_dxs());
    printf("number of hard particles: %d\n", get_hard());
    printf("number of emissions: %d\n", get_emissions());
    printf("\n");
    printf("particles:\n");
    for (int i = 0; i < get_size(); i++) {
      particle particle = get_particle(i);
      printf("  particle %d:\n", i + 1);
      printf("    pid: %d\n", particle.get_pid());
      printf("    mom: (%f, %f, %f, %f)\n", particle.get_mom()[0],
             particle.get_mom()[1], particle.get_mom()[2],
             particle.get_mom()[3]);
      printf("    [col, acol]: [%d, %d]\n", particle.get_col(),
             particle.get_acol());
      printf("    eta: %f\n", particle.get_eta());
    }
    printf("\n");
  }

  double gen_random() {
    /**
     * @brief generate a random number based on the event seed
     *
     * @return the random number
     */

    return generate_lcg(seed);
  }

  bool validate() {
    /**
     * @brief validate the event - check momentum and colour conservation
     *
     * @return bool: the validity of the event
     */

    vec4 psum = vec4();

    // n colours = n particles - 1
    int csum[max_particles - 1] = {0};

    // loop over particles
    for (int i = 0; i < get_size(); i++) {
      particle p = get_particle(i);

      // flip the momentum if it's an initial state particle
      vec4 pmom = p.is_initial() ? -p.get_mom() : p.get_mom();
      int pcol = p.get_col();
      int p_anti_col = p.get_acol();

      psum = psum + pmom;

      if (pcol > 0) {
        csum[pcol] += 1;
      }

      if (p_anti_col > 0) {
        csum[p_anti_col] -= 1;
      }
    }

    // Momentum Conservation: limit imbalance to 1e-7 GeV (10000 eV)
    double diff = 1e-5;
    bool pcheck = (fabs(psum[0]) < diff && fabs(psum[1]) < diff &&
                   fabs(psum[2]) < diff && fabs(psum[3]) < diff);

    // Colour Conservation
    bool ccheck = true;
    int broken_col = 0;
    for (int i = 0; i < max_particles - 1; i++) {
      if (csum[i] != 0) {
        ccheck = false;
        broken_col = i;
        break;
      }
    }

    // Output for failed validation
    if (!(pcheck && ccheck)) {
      std::cout << "event validation failed." << std::endl;

      if (!pcheck) {
        std::cout << "Momentum Imbalence: " << psum << std::endl;
      }

      if (!ccheck) {
        std::cout << "  Imbalence of Colour: " << broken_col << std::endl;
      }

      // print_info();
    }

    // Return validity
    return pcheck && ccheck;
  }
};

#endif  // event_h_