#ifndef chargedmult_h_
#define chargedmult_h_

#include "event.h"

// ---------------------------------------------------------------------------
// Electric charge from the digits of a PDG id (PDG Review, "Monte Carlo
// particle numbering scheme"):
//
//   meson : +/- nq2 nq3 nJ           (nq2 >= nq3, nJ = 2J+1)
//   baryon:     nq1 nq2 nq3 nJ       (nq1 >= nq2 >= nq3, nJ = 2J+1)

// ---------------------------------------------------------------------------
// Quark electric charge, in units of e/3 (up-type = +2, down-type = -1)

inline int qcharge3(int flav) {
  if (flav <= 0 || flav > 6) return 0;
  return (flav % 2 == 0) ? 2 : -1;
}

// ---------------------------------------------------------------------------
// Electric charge from PDG pid (antiparticle sign respected)

inline int particle_to_charge(int pid) {
  int a = abs(pid);
  int sign = (pid > 0) ? 1 : -1;

  // Leptons and gauge bosons
  switch (a) {
    case 11:
    case 13:
    case 15: return -sign;  // e, mu, tau
    case 12:
    case 14:
    case 16: return 0;      // neutrinos
    case 24: return sign;   // W+/-
    case 21:
    case 22:
    case 23:
    case 25: return 0;      // gluon, photon, Z, Higgs
  }

  int nq1 = (a / 1000) % 10;
  int nq2 = (a / 100) % 10;
  int nq3 = (a / 10) % 10;

  int c3;
  if (nq1 == 0) {
    // Meson: the heavier digit nq2 is a quark if it is up-type, else an
    // antiquark; nq3 is the opposite.  Charge = q(nq2) - q(nq3), oriented.
    int orient = (nq2 % 2 == 0) ? 1 : -1;
    c3 = orient * (qcharge3(nq2) - qcharge3(nq3));
  } else {
    // Baryon: sum of the three quark charges
    c3 = qcharge3(nq1) + qcharge3(nq2) + qcharge3(nq3);
  }

  return sign * c3 / 3;
}

// ---------------------------------------------------------------------------
// The charged multiplicity

void calculate_chargedmult(const event& ev, double* results);

#endif  // chargedmult_h_
