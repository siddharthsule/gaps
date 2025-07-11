#ifndef particle_cuh_
#define particle_cuh_

// particles have vec4 momentum, vec4 #includes base
#include "vec4.cuh"

class particle {
  /**
   * @class particle
   * @brief The particle object.
   *
   * This class is used to store the data for a particle. It contains the
   * particle id, the particle momentum, the particle colour, the particle
   * anti-colour, and the particle pseudorapidity.
   */

 public:
  // ---------------------------------------------------------------------------
  // member variables

  int pid;
  vec4 mom;
  int col;
  int acol;
  double eta;

  // ---------------------------------------------------------------------------
  // constructor

  __device__ particle(int pid = 0, vec4 momentum = vec4(), int col = 0,
                      int acol = 0, double eta = -1.)
      : pid(pid), mom(momentum), col(col), acol(acol), eta(eta) {}

  // ---------------------------------------------------------------------------
  // getters and setters
  __host__ __device__ int get_pid() const {
    /**
     * @brief get the particle id
     */

    return pid;
  }

  __host__ __device__ vec4 get_mom() const {
    /**
     * @brief get the particle momentum
     */

    return mom;
  }

  __host__ __device__ int get_col() const {
    /**
     * @brief get the particle colour
     */

    return col;
  }

  __host__ __device__ int get_acol() const {
    /**
     * @brief get the particle anti-colour
     */

    return acol;
  }

  __host__ __device__ double get_eta() const {
    /**
     * @brief get the particle pseudorapidity
     */

    return eta;
  }

  // ---------------------------------------------------------------------------
  // setters

  __device__ void set_pid(int pid) {
    /**
     * @brief set the particle id
     *
     * @param pid the particle id
     */

    this->pid = pid;
  }

  __device__ void set_mom(vec4 mom) {
    /**
     * @brief set the particle momentum
     *
     * @param mom the particle momentum
     */

    this->mom = mom;
  }

  __device__ void set_col(int col) {
    /**
     * @brief set the particle colour
     *
     * @param col the particle colour
     */

    this->col = col;
  }

  __device__ void set_acol(int acol) {
    /**
     * @brief set the particle anti-colour
     *
     * @param acol the particle anti-colour
     */

    this->acol = acol;
  }

  __device__ void set_eta(double eta) {
    /**
     * @brief set the particle pseudorapidity
     *
     * @param eta the particle pseudorapidity
     */

    this->eta = eta;
  }

  // ---------------------------------------------------------------------------
  // member functions

  __device__ bool is_initial() const {
    /**
     * @brief check if the particle is an initial state particle. If the
     * particle has a negative pseudorapidity, it is not an initial state
     * particle.
     */

    return eta < 0. ? false : true;
  }

  __device__ bool is_parton() const {
    /**
     * @brief check if the particle is a parton. If the particle id is 21 or
     * the abs value of the particle id is between 1 and 5, the particle is a
     * parton.
     */

    return (pid == 21) || (1 <= abs(pid) && abs(pid) <= 5);
  }

  // if two particles are in a colour connected dipole
  __device__ bool is_color_connected(particle p) {
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

#endif  // particle_cuh_