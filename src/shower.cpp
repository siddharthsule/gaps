#include "shower.h"

#include <cmath>
#include <random>

double Pqq::Value(double z, double y) {
  return kCF * (2. / (1. - z * (1. - y)) - (1. + z));
}

double Pqq::Estimate(double z) { return kCF * 2. / (1. - z); }

double Pqq::Integral(double zm, double zp) {
  return kCF * 2. * std::log((1. - zm) / (1. - zp));
}

double Pqq::GenerateZ(double zm, double zp) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  return 1. + (zp - 1.) * std::pow((1. - zm) / (1. - zp), dis(gen));
}

double Pgg::Value(double z, double y) {
  return kCA / 2. * (2. / (1. - z * (1. - y)) - 2. + z * (1. - z));
}

double Pgg::Estimate(double z) { return kCA / (1. - z); }

double Pgg::Integral(double zm, double zp) {
  return kCA * std::log((1. - zm) / (1. - zp));
}

double Pgg::GenerateZ(double zm, double zp) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  return 1. + (zp - 1.) * std::pow((1. - zm) / (1. - zp), dis(gen));
}

double Pgq::Value(double z, double y) {
  return kTR / 2. * (1. - 2. * z * (1. - z));
}

double Pgq::Estimate(double z) { return kTR / 2.; }

double Pgq::Integral(double zm, double zp) { return kTR / 2. * (zp - zm); }

double Pgq::GenerateZ(double zm, double zp) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  return zm + (zp - zm) * dis(gen);
}

Shower::Shower(double t0, AlphaS as) : t0(t0), as(as) {
  asmax = as(t0);

  int flavs[3];

  // 16 Kernels:
  // 10 Pqq, 1 Pgg, 5 Pgq
  int i = 0;
  for (int fl : {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5}) {
    flavs[0] = flavs[1] = fl;
    flavs[2] = 21;
    kernels[i] = new Pqq(flavs);
    i++;
  }

  flavs[0] = flavs[1] = flavs[2] = 21;
  kernels[i] = new Pgg(flavs);
  i++;

  for (int fl : {1, 2, 3, 4, 5}) {
    flavs[0] = 21;
    flavs[1] = fl;
    flavs[2] = -fl;
    kernels[i] = new Pgq(flavs);
    i++;
  }

}

void Shower::MakeKinematics(Vec4* kinematics, const double z, const double y,
                            const double phi, const Vec4 pijt, const Vec4 pkt) {
  Vec4 Q = pijt + pkt;

  // Generating the Momentum (0, kt1, kt2, 0)
  double rkt = sqrt(Q.M2() * y * z * (1. - z));

  Vec4 kt1 = pijt.Cross(pkt);
  if (kt1.P() < 1.e-6) {
    Vec4 xaxis(0., 1., 0., 0.);
    kt1 = pijt.Cross(xaxis);
  }
  kt1 = kt1 * (rkt * cos(phi) / kt1.P());

  Vec4 kt2cms = Q.Boost(pijt);
  kt2cms = kt2cms.Cross(kt1);
  kt2cms = kt2cms * (rkt * sin(phi) / kt2cms.P());
  Vec4 kt2 = Q.BoostBack(kt2cms);

  // Conversion to {i, j, k} basis
  Vec4 pi = pijt * z + pkt * ((1. - z) * y) + kt1 + kt2;
  Vec4 pj = pijt * (1. - z) + pkt * (z * y) - kt1 - kt2;
  Vec4 pk = pkt * (1. - y);

  kinematics[0] = pi;
  kinematics[1] = pj;
  kinematics[2] = pk;
}

void Shower::MakeColours(int* coli, int* colj, const int flavs[3],
                         const int colij[2], const int colk[2]) {
  // Increase the public variable c by 1
  c++;

  if (flavs[0] != 21) {
    if (flavs[0] > 0) {
      coli[0] = c;
      coli[1] = 0;
      colj[0] = colij[0];
      colj[1] = c;
    } else {
      coli[0] = 0;
      coli[1] = c;
      colj[0] = c;
      colj[1] = colij[1];
    }
  } else {
    if (flavs[1] == 21) {
      if (colij[0] == colk[1]) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        if (colij[1] == colk[0] && dis(gen) > 0.5) {
          coli[0] = colij[0];
          coli[1] = c;
          colj[0] = c;
          colj[1] = colij[1];
        } else {
          coli[0] = c;
          coli[1] = colij[1];
          colj[0] = colij[0];
          colj[1] = c;
        }
      } else {
        coli[0] = colij[0];
        coli[1] = c;
        colj[0] = c;
        colj[1] = colij[1];
      }
    } else {
      if (flavs[1] > 0) {
        coli[0] = colij[0];
        coli[1] = 0;
        colj[0] = 0;
        colj[1] = colij[1];
      } else {
        coli[0] = 0;
        coli[1] = colij[1];
        colj[0] = colij[0];
        colj[1] = 0;
      }
    }
  }
}

void Shower::GeneratePoint(Event& ev) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  // Winner parameters
  int win_ij = 0, win_k = 0;
  Kernel* win_sf;
  double win_m2 = 0.0, win_zp = 0.0;

  while (t > t0) {
    double t_temp = t0;

    for (size_t ij = 2; ij < ev.GetSize(); ij++) {
      for (size_t k = 2; k < ev.GetSize(); k++) {
        if (k == ij) {
          continue;
        }
        if (!ev.GetParton(ij).IsColorConnected(ev.GetParton(k))) {
          continue;
        }

        for (auto& sf : kernels) {
          if (sf->flavs[0] != ev.GetParton(ij).GetPid()) {
            continue;
          }

          double m2 = (ev.GetParton(ij).GetMom() + ev.GetParton(k).GetMom()).M2();
          if (m2 < 4.0 * t0) {
            continue;
          }
          double zp = 0.5 * (1.0 + std::sqrt(1.0 - 4.0 * t0 / m2));

          double g = asmax / (2.0 * M_PI) * sf->Integral(1.0 - zp, zp);
          double tt = t * std::pow(dis(gen), 1.0 / g);

          if (tt > t_temp) {
            t_temp = tt;
            win_ij = ij;
            win_k = k;
            win_sf = sf;
            win_m2 = m2;
            win_zp = zp;
          }
        }
      }
    }

    t = t_temp;

    if (t > t0) {
      double z = win_sf->GenerateZ(1.0 - win_zp, win_zp);
      double y = t / win_m2 / z / (1.0 - z);

      if (y < 1.0) {
        double f = (1.0 - y) * as(t) * win_sf->Value(z, y);
        double g = asmax * win_sf->Estimate(z);

        if (f / g > dis(gen)) {
          double phi = 2.0 * M_PI * dis(gen);

          Vec4 moms[3] = {Vec4(), Vec4(), Vec4()};
          MakeKinematics(moms, z, y, phi, ev.GetParton(win_ij).GetMom(),
                         ev.GetParton(win_k).GetMom());

          int colij[2] = {ev.GetParton(win_ij).GetCol(),
                          ev.GetParton(win_ij).GetAntiCol()};

          int colk[2] = {ev.GetParton(win_k).GetCol(),
                         ev.GetParton(win_k).GetAntiCol()};

          int coli[2] = {0, 0};
          int colj[2] = {0, 0};
          MakeColours(coli, colj, win_sf->flavs, colij, colk);

          ev.SetPartonPid(win_ij, win_sf->flavs[1]);
          ev.SetPartonMom(win_ij, moms[0]);
          ev.SetPartonCol(win_ij, coli[0]);
          ev.SetPartonAntiCol(win_ij, coli[1]);

          ev.SetPartonPid(ev.GetSize(), win_sf->flavs[2]);
          ev.SetPartonMom(ev.GetSize(), moms[1]);
          ev.SetPartonCol(ev.GetSize(), colj[0]);
          ev.SetPartonAntiCol(ev.GetSize(), colj[1]);

          ev.SetPartonMom(win_k, moms[2]);

          ev.IncrementEmissions();

          return;
        }
      }
    }
  }
}

void Shower::Run(Event& ev, double t_in) {
  t = t_in;
  c = 1;

  while (t > t0) {
    GeneratePoint(ev);
  }
}
