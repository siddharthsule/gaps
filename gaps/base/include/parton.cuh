#ifndef parton_cuh_
#define parton_cuh_

// partons have vec4 momentum, vec4 #includes base
#include "vec4.cuh"

/**
 * the parton class
 * ----------------

 * this file contains the parton object, which has attributes id, momentum and
 * colour. for now we use it for electrons too.
 */

class parton {
 public:
  // constructor
  __device__ parton(int pid = 0, vec4 momentum = vec4(), int col = 0,
                    int anticol = 0)
      : pid(pid), mom(momentum), col(col), anticol(anticol) {}

  // getters and setters
  __device__ int get_pid() const { return pid; }
  __device__ vec4 get_mom() const { return mom; }
  __device__ int get_col() const { return col; }
  __device__ int get_anti_col() const { return anticol; }

  __device__ void set_pid(int pid) { this->pid = pid; }
  __device__ void set_mom(vec4 mom) { this->mom = mom; }
  __device__ void set_col(int col) { this->col = col; }
  __device__ void set_anti_col(int anticol) { this->anticol = anticol; }

  // if two partons are in a colour connected dipole
  __device__ bool is_color_connected(parton p) {
    return (col > 0 && col == p.anticol) || (anticol > 0 && anticol == p.col);
  }

 private:
  int pid;
  vec4 mom;
  int col;
  int anticol;
};

#endif  // parton_cuh_