#ifndef database_cuh_
#define database_cuh_

#include "base.cuh"

// ---------------------------------------------------------------------------
// Hadron masses [GeV], lifetimes [mm] and decay modes, keyed by
// signed pid
//
// GENERATED FILE - do not edit by hand.

// ctau of a final state, and of a state with no PDG width (K0,
// K0bar); neither is a lifetime to compare against ctau_max.
inline constexpr double ctau_stable = 1e30;
inline constexpr double ctau_no_width = -1.;

// The most children a mode may name
inline constexpr int max_decay_products = 10;

// Struct to hold hadron information
struct hadron_row {
  int pid;         // signed, as the PDG numbers it
  double mass;     // GeV
  double ctau;     // mm, or one of the sentinels above
  int mode_begin;  // first row of its modes in decay_database
  int n_modes;     // how many, zero for a state that decays nowhere
};

// One child: K0/K0bar mixing into K_S/K_L, with all the momentum
struct decay_row {
  double br;                         // branching fraction, unnormalised
  int n_children;                    // 1 to max_decay_products
  int children[max_decay_products];  // their pids, the rest unused
};

inline constexpr int n_hadron_rows = 265;
inline constexpr int n_decay_rows = 6383;

// ---------------------------------------------------------------------------
// The tables, defined in database.cu

extern __device__ const hadron_row hadron_database[n_hadron_rows];
extern __device__ const decay_row decay_database[n_decay_rows];

// ---------------------------------------------------------------------------
// Row of a state in hadron_database, or < 0 where the table does
// not hold it.

__device__ inline int hadron_index(int pid) {
  int lo = 0;
  int hi = n_hadron_rows - 1;

  while (lo <= hi) {
    int mid = (lo + hi) / 2;
    if (hadron_database[mid].pid == pid) return mid;
    if (hadron_database[mid].pid < pid) {
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }

  return -1;
}

// ---------------------------------------------------------------------------
// Mass of a state [GeV], or < 0 where the table does not hold
// it.

__device__ inline double hadron_mass(int pid) {
  int i = hadron_index(pid);
  return (i < 0) ? -1.0 : hadron_database[i].mass;
}

#endif  // database_cuh_
