#ifndef particle_h_
#define particle_h_

#include "vec4.h"

class particle {
  /**
   * @class particle
   * @brief The particle object.
   *
   * This class is used to store the data for a particle. It contains the
   * particle id, the particle momentum, the particle colour, the particle
   * anti-colour, and the particle pseudorapidity.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  int pid;
  vec4 mom;
  int col;
  int acol;
  double eta;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  particle(int pid = 0, vec4 momentum = vec4(), int col = 0, int acol = 0,
           double eta = -1.)
      : pid(pid), mom(momentum), col(col), acol(acol), eta(eta) {}

  // ---------------------------------------------------------------------------
  // getters

  int get_pid() const {
    /**
     * @brief get the particle id
     */

    return pid;
  }

  vec4 get_mom() const {
    /**
     * @brief get the particle momentum
     */

    return mom;
  }

  int get_col() const {
    /**
     * @brief get the particle colour
     */

    return col;
  }

  int get_acol() const {
    /**
     * @brief get the particle anti-colour
     */

    return acol;
  }

  double get_eta() const {
    /**
     * @brief get the particle pseudorapidity
     */

    return eta;
  }

  // ---------------------------------------------------------------------------
  // setters

  void set_pid(int pid) {
    /**
     * @brief set the particle id
     *
     * @param pid the particle id
     */

    this->pid = pid;
  }

  void set_mom(vec4 mom) {
    /**
     * @brief set the particle momentum
     *
     * @param mom the particle momentum
     */

    this->mom = mom;
  }

  void set_col(int col) {
    /**
     * @brief set the particle colour
     *
     * @param col the particle colour
     */

    this->col = col;
  }

  void set_acol(int acol) {
    /**
     * @brief set the particle anti-colour
     *
     * @param acol the particle anti-colour
     */

    this->acol = acol;
  }

  void set_eta(double eta) {
    /**
     * @brief set the particle pseudorapidity
     *
     * @param eta the particle pseudorapidity
     */

    this->eta = eta;
  }

  // ---------------------------------------------------------------------------
  // member functions

  bool is_initial() const {
    /**
     * @brief check if the particle is an initial state particle. If the
     * particle has a negative pseudorapidity, it is not an initial state
     * particle.
     */

    return eta < 0. ? false : true;
  }

  bool is_parton() const {
    /**
     * @brief check if the particle is a parton. If the particle id is 21 or
     * the abs value of the particle id is between 1 and 5, the particle is a
     * parton.
     */

    return (pid == 21) || (1 <= abs(pid) && abs(pid) <= 5);
  }

  bool is_color_connected(particle p) {
    /**
     * @brief check if the particle is in a colour connected dipole with another
     * particle. If the particle colour is the same as the other particle
     * anti-colour or the particle anti-colour is the same as the other particle
     * colour, the particles are in a colour connected dipole.
     *
     * @param p the other particle
     */

    return (col > 0 && col == p.acol) || (acol > 0 && acol == p.col);
  }
};

#endif  // particle_h_
