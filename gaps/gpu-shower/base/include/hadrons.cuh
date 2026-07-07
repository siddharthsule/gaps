#ifndef hadrons_cuh_
#define hadrons_cuh_

#include "base.cuh"

// ---------------------------------------------------------------------------
// PDG Monte Carlo particle numbering — computed, not tabulated
//
// Rather than scanning a table of quark-content -> hadron rows, we build the
// PDG id straight from the constituent flavours using the numbering rules
// (PDG Review, "Monte Carlo particle numbering scheme"):
//
//   meson : +/- nq2 nq3 nJ           (nq2 >= nq3, nJ = 2J+1)
//   baryon:     nq1 nq2 nq3 nJ       (nq1 >= nq2 >= nq3, nJ = 2J+1)
//
// Only the *mass* is genuinely data (a physical, partly model-tuned constant),
// so that is the one thing kept in a lookup — keyed by |pid|, one row per
// distinct hadron.  Charge is also derived from the digits.
//
// Computing the ids removes the __constant__ lookup tables entirely: every
// thread evaluates the same branch-light arithmetic, so there is no table to
// broadcast and no separable-compilation multiply-definition to worry about.
//
// Sign conventions here follow strict PDG: e.g. (s, dbar) is Kbar0 = -311, the
// mirror of (d, sbar) = K0 = +311.
//
// Flavour codes: d=1, u=2, s=3, c=4, b=5.

// ---------------------------------------------------------------------------
// Quark electric charge, in units of e/3 (up-type = +2, down-type = -1)

__device__ inline int qcharge3(int flav) {
  if (flav <= 0 || flav > 6) return 0;
  return (flav % 2 == 0) ? 2 : -1;
}

// ---------------------------------------------------------------------------
// Diquark identification
//
// A diquark id has the form nq1 nq2 0 nJ with nq1,nq2 >= 1 and nJ in {1,3}
// (scalar / vector).  This is exactly the PDG structure, so we test the digits
// instead of listing every diquark explicitly.

__device__ inline bool is_diquark(int pid) {
  int a = abs(pid);
  if (a < 1000 || a >= 10000) return false;
  int nq1 = (a / 1000) % 10;
  int nq2 = (a / 100) % 10;
  int nq3 = (a / 10) % 10;
  int nj = a % 10;
  return nq3 == 0 && nq1 >= 1 && nq2 >= 1 && (nj == 1 || nj == 3);
}

// ---------------------------------------------------------------------------
// Hadron mass [GeV], keyed by |pid| (one row per distinct hadron)
//
// Physical values, except phi (333) whose mass is tuned to 0.890 for low-mass
// s sbar clusters in the hadron+gamma fallback.  Returns < 0 for any hadron
// this model does not produce (e.g. charm/bottom baryons), which is what gates
// the lookups below.

__device__ inline double hadron_mass(int pid) {
  switch (abs(pid)) {
    // light / strange mesons
    case 111: return 0.135;   // pi0
    case 211: return 0.139;   // pi+/-
    case 311: return 0.498;   // K0 / Kbar0
    case 321: return 0.494;   // K+/-
    case 333: return 0.890;   // phi (tuned)
    // charm / bottom mesons
    case 411: return 1.870;   // D+/-
    case 421: return 1.865;   // D0 / D0bar
    case 431: return 1.968;   // Ds+/-
    case 443: return 3.097;   // J/psi
    case 511: return 5.280;   // B0 / B0bar
    case 521: return 5.279;   // B+/-
    case 531: return 5.370;   // Bs0 / Bs0bar
    case 541: return 6.275;   // Bc+/-
    case 553: return 9.460;   // Upsilon
    // baryons
    case 2212: return 0.9383;  // proton
    case 2112: return 0.9396;  // neutron
    case 3122: return 1.1157;  // Lambda0
    case 3222: return 1.1894;  // Sigma+
    case 3212: return 1.1925;  // Sigma0
    case 3112: return 1.1974;  // Sigma-
    case 3322: return 1.3149;  // Xi0
    case 3312: return 1.3217;  // Xi-
    case 3334: return 1.6725;  // Omega-
    default: return -1.0;
  }
}

// ---------------------------------------------------------------------------
// Meson pid from a quark / antiquark pair
//
// q_pid > 0 (quark), qbar_pid < 0 (antiquark).  Returns the signed PDG id, or
// 0 if the two ids are not a quark-antiquark pair among d,u,s,c,b.

__device__ inline int meson_pid(int q_pid, int qbar_pid) {
  if (q_pid < 0 || qbar_pid > 0) return 0;  // must be (quark, antiquark)
  int a = q_pid, b = -qbar_pid;             // flavours 1..5
  if (a < 1 || a > 5 || b < 1 || b > 5) return 0;

  // Flavour-diagonal (quarkonium): a model choice of the ground state kept —
  // pseudoscalar pi0 for light, vector for s/c/b (phi, J/psi, Upsilon).
  if (a == b) {
    switch (a) {
      case 1:
      case 2: return 111;  // pi0
      case 3: return 333;  // phi   (vector, tuned mass)
      case 4: return 443;  // J/psi
      case 5: return 553;  // Upsilon
    }
  }

  int fh = (a > b) ? a : b;  // heavier flavour -> nq2
  int fl = (a < b) ? a : b;  // lighter flavour -> nq3
  int code = 100 * fh + 10 * fl + 1;  // ground-state pseudoscalar (2J+1 = 1)

  // Sign: a meson is a particle (+) if its heavier quark is an up-type quark
  // or a down-type antiquark.  h is the signed id of the heavier flavour.
  int h = (a >= b) ? q_pid : qbar_pid;
  int sign = (fh % 2 == 0) ? ((h > 0) ? 1 : -1) : ((h > 0) ? -1 : 1);
  return sign * code;
}

// ---------------------------------------------------------------------------
// Baryon pid from a quark + diquark pair
//
// q_pid is a (anti)quark, diq_pid a (anti)diquark of the same sign.  Returns
// the signed PDG id (matter baryon for quark > 0), or 0 on bad input.

__device__ inline int baryon_pid(int q_pid, int diq_pid) {
  int sign = (q_pid > 0) ? 1 : -1;
  int qf = abs(q_pid);
  int d = abs(diq_pid);
  int d1 = (d / 1000) % 10;
  int d2 = (d / 100) % 10;
  bool scalar_diq = (d % 10 == 1);  // spin-0 diquark
  if (qf < 1 || d1 < 1 || d2 < 1) return 0;

  // Sort the three flavours descending -> nq1 >= nq2 >= nq3
  int n1 = qf, n2 = d1, n3 = d2;
  if (n1 < n2) { int t = n1; n1 = n2; n2 = t; }
  if (n2 < n3) { int t = n2; n2 = n3; n3 = t; }
  if (n1 < n2) { int t = n1; n1 = n2; n2 = t; }

  // Ground-state spin: three identical flavours are forced to J=3/2 (Pauli,
  // e.g. Omega-, Delta), everything else is J=1/2.
  int nj = (n1 == n2 && n2 == n3) ? 4 : 2;

  // Lambda vs Sigma0 for three distinct flavours: the lighter (Lambda-type)
  // state has the two light quarks in a spin-0 diquark, with nq2/nq3 swapped.
  // Here that means the free quark is the heaviest and the diquark is scalar.
  bool all_distinct = (n1 != n2 && n2 != n3);
  bool lambda_type = all_distinct && scalar_diq && (qf == n1);

  int base = lambda_type ? (1000 * n1 + 100 * n3 + 10 * n2 + 2)
                         : (1000 * n1 + 100 * n2 + 10 * n3 + nj);
  return sign * base;
}

// ---------------------------------------------------------------------------
// Meson lookup — sets out_pid/out_mass, returns false if unsupported.

__device__ inline bool which_meson(int q_pid, int qbar_pid, int& out_pid,
                                   double& out_mass) {
  int pid = meson_pid(q_pid, qbar_pid);
  if (pid == 0) return false;
  double m = hadron_mass(pid);
  if (m < 0) return false;
  out_pid = pid;
  out_mass = m;
  return true;
}

// ---------------------------------------------------------------------------
// Baryon lookup
//
// One of (pid1, pid2) is a (anti)diquark, the other is a (anti)quark.

__device__ inline bool which_baryon(int pid1, int pid2, int& out_pid,
                                    double& out_mass) {
  int q_pid, diq_pid;
  if (is_diquark(pid1)) {
    diq_pid = pid1;
    q_pid = pid2;
  } else {
    q_pid = pid1;
    diq_pid = pid2;
  }

  int pid = baryon_pid(q_pid, diq_pid);
  if (pid == 0) return false;
  double m = hadron_mass(pid);
  if (m < 0) return false;
  out_pid = pid;
  out_mass = m;
  return true;
}

// ---------------------------------------------------------------------------
// Unified hadron lookup

__device__ inline bool which_hadron(int pid1, int pid2, int& out_pid,
                                    double& out_mass) {
  if (is_diquark(pid1) || is_diquark(pid2))
    return which_baryon(pid1, pid2, out_pid, out_mass);
  return which_meson(pid1, pid2, out_pid, out_mass);
}

// ---------------------------------------------------------------------------
// Electric charge from PDG pid (antiparticle sign respected)
//
// Leptons and gauge bosons are handled explicitly; hadrons are summed from
// their quark-content digits (meson = quark - antiquark, baryon = 3 quarks).

__device__ inline int particle_to_charge(int pid) {
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

#endif  // hadrons_cuh_
