#ifndef PARTON_H_
#define PARTON_H_

// Partons have Vec4 Momentum, Vec4 #includes Base
#include "vec4.h"

class Parton {
 public:
  // Constructor
  // Used at Multiple Points [HOST + DEVICE]
  Parton(int pid = 0, Vec4 momentum = Vec4(), int col = 0, int anticol = 0,
         bool initial = false, bool parton = true)
      : pid(pid),
        mom(momentum),
        col(col),
        anticol(anticol),
        initial(initial),
        parton(parton) {}

  // Getters and Setters
  int GetPid() const { return pid; }
  Vec4 GetMom() const { return mom; }
  int GetCol() const { return col; }
  int GetAntiCol() const { return anticol; }
  bool IsInitial() const { return initial; }
  bool IsParton() const { return parton; }

  void SetPid(int pid) { this->pid = pid; }
  void SetMom(Vec4 mom) { this->mom = mom; }
  void SetCol(int col) { this->col = col; }
  void SetAntiCol(int anticol) { this->anticol = anticol; }
  void SetInitial(bool initial) { this->initial = initial; }
  void SetParton(bool parton) { this->parton = parton; }

  // Boolean - If two partons are in a Colour Connected Dipole
  bool IsColorConnected(Parton p) {
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

#endif  // PARTON_H_
