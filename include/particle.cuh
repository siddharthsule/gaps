#ifndef PARTICLE_CUH_
#define PARTICLE_CUH_

#include <cuda_runtime.h>

// Particles have Vec4 Momentum
#include "vec4.cuh"

class Particle {
 public:
  // Constructor
  // Used at Multiple Points [HOST + DEVICE]
  __host__ __device__ Particle(int pid = 0, Vec4 momentum = Vec4(), int col = 0,
                               int anticol = 0)
      : pid(pid), mom(momentum), col(col), anticol(anticol) {}

  // Getters and Setters
  // Not used in ME [HOST]
  int GetPid() const { return pid; }
  Vec4 GetMom() const { return mom; }
  int GetCol() const { return col; }
  int GetAntiCol() const { return anticol; }

  void SetPid(int pid) { this->pid = pid; }
  void SetMom(Vec4 mom) { this->mom = mom; }
  void SetCol(int col) { this->col = col; }
  void SetAntiCol(int anticol) { this->anticol = anticol; }

  // If two partons are in a Colour Connected Pair
  // For Parton Shower [HOST]
  // Boolean - If two partons are in a Colour Connected Pair
  bool IsColorConnected(Particle p) {
    return (col > 0 && col == p.anticol) || (anticol > 0 && anticol == p.col);
  }

 private:
  int pid;
  Vec4 mom;
  int col;
  int anticol;
};

// Event Class
// Built to contain the partons and the dxs as one accesible object
// In future, can use to store thrust, log10y23 to parallelise those
class Event {
 private:
  // Very Temporary Solution - Allows 30 Emissions
  // Mainly done to avoid using std::vector for GPU Showering
  // But of course, this is not a good solution
  Particle partons[34];
  double dxs;
  const int nHard = 4;
  int nEmissions;

 public:
  // Constructor
  __host__ __device__ Event(Particle input_p[4], double dxs) : dxs(dxs) {
    for (int i = 0; i < 4; i++) {
      partons[i] = input_p[i];
    }
    nEmissions = 0;
  }

  // Getters
  // Not used in ME [HOST]
  Particle GetParton(int i) const { return partons[i]; }
  double GetDxs() const { return dxs; }

  // Get Total Number of Particles
  int GetSize() const { return nHard + nEmissions; }

  // Setters
  // Not used in ME [HOST]
  void SetPartonPid(int i, int pid) { partons[i].SetPid(pid); }
  void SetPartonMom(int i, Vec4 mom) { partons[i].SetMom(mom); }
  void SetPartonCol(int i, int col) { partons[i].SetCol(col); }
  void SetPartonAntiCol(int i, int anticol) { partons[i].SetAntiCol(anticol); }

  // Increment Number of Emissions
  void IncrementEmissions() { nEmissions++; }
};

#endif  // PARTICLE_CUH_