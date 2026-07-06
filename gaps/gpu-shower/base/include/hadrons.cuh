#ifndef hadrons_cuh_
#define hadrons_cuh_

#include "base.cuh"

// ---------------------------------------------------------------------------
// Diquark identification

__device__ inline bool is_diquark(int pid) {
  int a = abs(pid);
  return (a == 2101 || a == 2103 || a == 2203 || a == 1103 || a == 3101 ||
          a == 3103 || a == 3201 || a == 3203 || a == 3303);
}

// ---------------------------------------------------------------------------
// Hadron lookup tables
//
// MesonEntry  — q > 0 (quark pid), qbar < 0 (antiquark pid).
//   pid is the signed PDG id; charge is the electric charge for that pid.
//
// BaryonEntry — q = |quark pid|, diq = |diquark pid|.
//   pid is always the matter-baryon PDG id (> 0); charge is the charge of
//   the matter baryon.  The sign of the output pid is determined from the
//   sign of the input quark at lookup time.
//
// The tables live in __constant__ memory: all threads in a warp scan the
// same entry in lockstep, which the constant cache serves as a broadcast.
// `static` gives each translation unit its own private copy, as required
// under separable compilation (a shared-header __constant__ with external
// linkage would multiply-define at device-link time).

struct MesonEntry {
  int q, qbar;  // quark / antiquark pids
  int pid;      // signed PDG pid
  double mass;  // mass [GeV]
  int charge;   // electric charge for this pid
};

struct BaryonEntry {
  int q, diq;   // |quark pid|, |diquark pid|
  int pid;      // matter-baryon PDG pid (always > 0)
  double mass;  // mass [GeV]
  int charge;   // electric charge of the matter baryon
};

static __constant__ MesonEntry meson_table[] = {
    {2, -1, 211, 0.139, 1},    // pi+      (u dbar)
    {1, -2, -211, 0.139, -1},  // pi-      (d ubar)
    {1, -1, 111, 0.135, 0},    // pi0      (d dbar)
    {2, -2, 111, 0.135, 0},    // pi0      (u ubar)
    {2, -3, 321, 0.494, 1},    // K+       (u sbar)
    {3, -2, -321, 0.494, -1},  // K-       (s ubar)
    {1, -3, 311, 0.498, 0},    // K0       (d sbar)
    {3, -1, 311, 0.498, 0},    // K0       (s dbar)
    {3, -3, 333, 0.890, 0},    // phi      (s sbar), tuned for low-mass s sbar
                               // clusters in hadron+gamma fallback
    {4, -1, 411, 1.870, 1},    // D+       (c dbar)
    {1, -4, -411, 1.870, -1},  // D-       (d cbar)
    {4, -2, 421, 1.865, 0},    // D0       (c ubar)
    {2, -4, 421, 1.865, 0},    // D0       (u cbar)
    {4, -3, 431, 1.968, 1},    // Ds+      (c sbar)
    {3, -4, -431, 1.968, -1},  // Ds-      (s cbar)
    {4, -4, 443, 3.097, 0},    // J/psi    (c cbar)
    {2, -5, 521, 5.279, 1},    // B+       (u bbar)
    {5, -2, -521, 5.279, -1},  // B-       (b ubar)
    {1, -5, 511, 5.280, 0},    // B0       (d bbar)
    {5, -1, 511, 5.280, 0},    // B0       (b dbar)
    {3, -5, 531, 5.370, 0},    // Bs0      (s bbar)
    {5, -3, 531, 5.370, 0},    // Bs0      (b sbar)
    {4, -5, 541, 6.275, 1},    // Bc+      (c bbar)
    {5, -4, -541, 6.275, -1},  // Bc-      (b cbar)
    {5, -5, 553, 9.460, 0},    // Upsilon  (b bbar)
};

static __constant__ BaryonEntry baryon_table[] = {
    {2, 2101, 2212, 0.9383, 1},   // proton   (u ud)
    {2, 2103, 2212, 0.9383, 1},   // proton   (u ud)
    {1, 2203, 2212, 0.9383, 1},   // proton   (u ud)
    {1, 2101, 2112, 0.9396, 0},   // neutron  (d ud)
    {1, 2103, 2112, 0.9396, 0},   // neutron  (d ud)
    {2, 1103, 2112, 0.9396, 0},   // neutron  (d ud)
    {3, 2101, 3122, 1.1157, 0},   // Lambda0  (s ud)
    {2, 3103, 3222, 1.1894, 1},   // Sigma+   (u us)
    {3, 2203, 3222, 1.1894, 1},   // Sigma+   (u us)
    {2, 3203, 3212, 1.1925, 0},   // Sigma0   (u ds)
    {1, 3103, 3212, 1.1925, 0},   // Sigma0   (d us)
    {3, 2103, 3212, 1.1925, 0},   // Sigma0   (d us)
    {1, 3203, 3112, 1.1974, -1},  // Sigma-   (d ds)
    {3, 1103, 3112, 1.1974, -1},  // Sigma-   (d ds)
    {2, 3303, 3322, 1.3149, 0},   // Xi0      (u ss)
    {3, 3103, 3322, 1.3149, 0},   // Xi0      (u ss)
    {1, 3303, 3312, 1.3217, -1},  // Xi-      (d ss)
    {3, 3203, 3312, 1.3217, -1},  // Xi-      (d ss)
    {3, 3303, 3334, 1.6725, -1},  // Omega-   (s ss)
};

// ---------------------------------------------------------------------------
// Meson lookup
//
// q_pid > 0 (quark), qbar_pid < 0 (antiquark).
// Returns true and sets out_pid/out_mass on success, false otherwise.

__device__ inline bool which_meson(int q_pid, int qbar_pid, int& out_pid,
                                   double& out_mass) {
  for (const auto& e : meson_table) {
    if (e.q == q_pid && e.qbar == qbar_pid) {
      out_pid = e.pid;
      out_mass = e.mass;
      return true;
    }
  }
  return false;
}

// ---------------------------------------------------------------------------
// Baryon lookup
//
// One of (pid1, pid2) is a (anti)diquark, the other is a (anti)quark.
// Returns true and sets out_pid/out_mass on success, false otherwise.

__device__ inline bool which_baryon(int pid1, int pid2, int& out_pid,
                                    double& out_mass) {
  int q_pid, diq_pid;
  if (is_diquark(abs(pid1))) {
    diq_pid = pid1;
    q_pid = pid2;
  } else {
    q_pid = pid1;
    diq_pid = pid2;
  }

  int sign = (q_pid > 0) ? 1 : -1;
  int q = abs(q_pid);
  int diq = abs(diq_pid);

  for (const auto& e : baryon_table) {
    if (e.q == q && e.diq == diq) {
      out_pid = sign * e.pid;
      out_mass = e.mass;
      return true;
    }
  }
  return false;
}

// ---------------------------------------------------------------------------
// Unified hadron lookup

__device__ inline bool which_hadron(int pid1, int pid2, int& out_pid,
                                    double& out_mass) {
  if (is_diquark(abs(pid1)) || is_diquark(abs(pid2)))
    return which_baryon(pid1, pid2, out_pid, out_mass);
  return which_meson(pid1, pid2, out_pid, out_mass);
}

// ---------------------------------------------------------------------------
// Electric charge from PDG pid (antiparticle sign is respected)
//
// For mesons: uses the first table entry whose stored pid is positive and
// matches |pid|, giving the canonical charge; sign(pid) is then applied.
// For baryons: all table entries store the matter-baryon (positive) pid;
// sign(pid) is applied the same way.

__device__ inline int particle_to_charge(int pid) {
  int a = pid > 0 ? pid : -pid;

  for (const auto& e : meson_table) {
    if (e.pid > 0 && e.pid == a) return pid > 0 ? e.charge : -e.charge;
  }

  for (const auto& e : baryon_table) {
    if (e.pid == a) return pid > 0 ? e.charge : -e.charge;
  }

  // Leptons and gauge bosons
  const int lg_pids[] = {11, 13, 15, 22, 23, 24};
  const int lg_charges[] = {-1, -1, -1, 0, 0, 1};
  for (int i = 0; i < 6; ++i) {
    if (lg_pids[i] == a) return pid > 0 ? lg_charges[i] : -lg_charges[i];
  }

  return 0;
}

#endif  // hadrons_cuh_
