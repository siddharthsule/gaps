#ifndef SHOWER_H_
#define SHOWER_H_

#include <vector>

#include "particle.cuh"
#include "qcd.h"
#include "vec4.cuh"

class Kernel {
public:

  explicit Kernel(int fl[3]) {
    for (int i = 0; i < 3; i++) flavs[i] = fl[i];
  }

  virtual double Value(double z, double y) = 0;
  virtual double Estimate(double z) = 0;
  virtual double Integral(double zm, double zp) = 0;
  virtual double GenerateZ(double zm, double zp) = 0;

public:
  int flavs[3];
};

class Pqq : public Kernel {
 public:
  explicit Pqq(int fl[3]) : Kernel(fl){};

  double Value(double z, double y) override;
  double Estimate(double z) override;
  double Integral(double zm, double zp) override;
  double GenerateZ(double zm, double zp) override;
};

class Pgg : public Kernel {
 public:
  explicit Pgg(int fl[3]) : Kernel(fl){};

  double Value(double z, double y) override;
  double Estimate(double z) override;
  double Integral(double zm, double zp) override;
  double GenerateZ(double zm, double zp) override;
};

class Pgq : public Kernel {
 public:
  explicit Pgq(int fl[3]) : Kernel(fl){};

  double Value(double z, double y) override;
  double Estimate(double z) override;
  double Integral(double zm, double zp) override;
  double GenerateZ(double zm, double zp) override;
};

class Shower {
 public:
  explicit Shower(double t0, AlphaS as);

  void MakeKinematics(Vec4* kinematics, const double z, const double y,
                      const double phi, const Vec4 pijt, const Vec4 pkt);
  void MakeColours(int* coli, int* colj, const int flavs[3], const int colij[2],
                   const int colk[2]);
  void GeneratePoint(Event& ev);
  void Run(Event& ev, double t);

 private:
  int c;
  double t, t0;
  AlphaS as;
  double asmax;
  int nEmissions;
  Kernel* kernels[16];
};

#endif  // SHOWER_H_