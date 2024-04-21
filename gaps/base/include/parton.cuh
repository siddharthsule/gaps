#ifndef PARTON_CUH_
#define PARTON_CUH_

// Partons have Vec4 Momentum, Vec4 #includes Base
#include "vec4.cuh"

/**
 * The Parton Class
 * ----------------

 * This file contains the Parton Object, which has attributes ID, momentum and
 * colour. For now we use it for Electrons too.
 */

class Parton {
 public:
  // Constructor
  __device__ Parton(int pid = 0, Vec4 momentum = Vec4(), int col = 0,
                    int anticol = 0, bool initial = false, bool parton = true)
      : pid(pid),
        mom(momentum),
        col(col),
        anticol(anticol),
        initial(initial),
        parton(parton) {}

  // Getters and Setters
  __device__ int GetPid() const { return pid; }
  __device__ Vec4 GetMom() const { return mom; }
  __device__ int GetCol() const { return col; }
  __device__ int GetAntiCol() const { return anticol; }
  // Prep for Initial State (v2...)
  __device__ bool IsInitial() const { return initial; }
  __device__ bool IsParton() const { return parton; }

  __device__ void SetPid(int pid) { this->pid = pid; }
  __device__ void SetMom(Vec4 mom) { this->mom = mom; }
  __device__ void SetCol(int col) { this->col = col; }
  __device__ void SetAntiCol(int anticol) { this->anticol = anticol; }
  __device__ void SetInitial(bool initial) { this->initial = initial; }
  __device__ void SetParton(bool parton) { this->parton = parton; }

  // If two partons are in a Colour Connected Dipole
  __device__ bool IsColorConnected(Parton p) {
    return (col > 0 && col == p.anticol) || (anticol > 0 && anticol == p.col);
  }

 private:
  int pid;
  Vec4 mom;
  int col;
  int anticol;
  bool initial;
  bool parton;
};

#endif  // PARTON_CUH_