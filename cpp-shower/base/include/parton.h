#ifndef parton_h_
#define parton_h_

// partons have vec4 momentum, vec4 #includes base
#include "vec4.h"

class parton {
 public:
  // constructor
  parton(int pid = 0, vec4 momentum = vec4(), int col = 0, int anticol = 0)
      : pid(pid), mom(momentum), col(col), anticol(anticol) {}

  // getters and setters
  int get_pid() const { return pid; }
  vec4 get_mom() const { return mom; }
  int get_col() const { return col; }
  int get_anti_col() const { return anticol; }

  void set_pid(int pid) { this->pid = pid; }
  void set_mom(vec4 mom) { this->mom = mom; }
  void set_col(int col) { this->col = col; }
  void set_anti_col(int anticol) { this->anticol = anticol; }

  // boolean - if two partons are in a colour connected dipole
  bool is_color_connected(parton p) {
    return (col > 0 && col == p.anticol) || (anticol > 0 && anticol == p.col);
  }

 private:
  int pid;
  vec4 mom;
  int col;
  int anticol;
};

#endif  // parton_h_
