#ifndef VEC4_CUH_
#define VEC4_CUH_

// Base Class, with all the important definitions
#include "base.cuh"

/**
 * Four Momenta
 * ------------
 *
 * This file contains the definition of the Four Momenta class, which is used to
 * represent the four-momentum of particles in the event. It is a simple class
 * with the four attributes (E, px, py, pz) and some basic operations that can
 * be performed with them.
 *
 *
 * Why is this Header Only?
 * ------------------------
 * 
 * When you declare a __device__ function in a .cuh (header) file and define it 
 * in a .cu file, that function's definition is only available to the .cu file 
 * in which it is defined. This is because __device__ functions are compiled by 
 * nvcc into the device code, and unlike host functions, they do not have 
 * external linkage that allows them to be seen or linked across different .cu 
 * files after individual compilation.
 */

class Vec4 {
 private:
  double E, px, py, pz;

 public:
  // Constructor - Define key attributes Energy and Momentum
  // Used in ME and out [HOST + DEVICE]
  __device__ Vec4(double E = 0., double px = 0., double py = 0., double pz = 0.)
      : E(E), px(px), py(py), pz(pz) {}

  // Get Method to Obtain Attribute Value
  __device__ double operator[](int i) const {
    switch (i) {
      case 0:
        return E;
      case 1:
        return px;
      case 2:
        return py;
      case 3:
        return pz;
      default:
        // CUDA does not support exceptions, so we just return 0
        return 0;
    }
  }

  // Print a Column Vector with the attributes
  friend std::ostream& operator<<(std::ostream& os, const Vec4& v) {
    os << "(" << v.E << "," << v.px << "," << v.py << "," << v.pz << ")";
    return os;
  }

  // Simple Mathematics with Four vectors
  __device__ Vec4 operator+(const Vec4& v) const {
    return Vec4(E + v.E, px + v.px, py + v.py, pz + v.pz);
  }

  __device__ Vec4 operator-() const { return Vec4(-E, -px, -py, -pz); }

  __device__ Vec4 operator-(const Vec4& v) const {
    return Vec4(E - v.E, px - v.px, py - v.py, pz - v.pz);
  }

  // Multiplication (and Dot Product)
  __device__ double operator*(const Vec4& v) const {
    return E * v.E - px * v.px - py * v.py - pz * v.pz;
  }

  __device__ Vec4 operator*(double v) const {
    return Vec4(E * v, px * v, py * v, pz * v);
  }

  // Division
  __device__ Vec4 operator/(double v) const {
    return Vec4(E / v, px / v, py / v, pz / v);
  }

  // Magnitude of the Vector
  __device__ double M2() const { return (*this) * (*this); }

  __device__ double M() const {
    double m2 = M2();
    return m2 > 0 ? sqrt(m2) : 0;
  }

  // 3 Momenta
  __device__ double P2() const { return px * px + py * py + pz * pz; }

  __device__ double P() const {
    double p2 = P2();
    return p2 > 0 ? sqrt(p2) : 0;
  }

  // Transverse Momenta
  __device__ double PT2() const { return px * px + py * py; }

  __device__ double PT() const {
    double pt2 = PT2();
    return pt2 > 0 ? sqrt(pt2) : 0;
  }

  // Angles
  __device__ double Theta() const {
    double p = P();
    return p != 0 ? acos(pz / p) : 0;
  }

  __device__ double Phi() const {
    if (px == 0 && py == 0) {
      return 0.0;
    } else {
      return atan2(py, px);
    }
  }

  __device__ double Rapidity() const {
    double denominator = (E - pz);
    return denominator != 0 ? 0.5 * log((E + pz) / denominator) : 0;
  }

  __device__ double Eta() const {
    double theta = Theta();
    return -log(tan(theta / 2.));
  }

  // Three Momenta Dot and Cross Product
  __device__ double Dot(const Vec4& v) const {
    return px * v.px + py * v.py + pz * v.pz;
  }

  __device__ Vec4 Cross(const Vec4& v) const {
    return Vec4(0.0, py * v.pz - pz * v.py, pz * v.px - px * v.pz,
                px * v.py - py * v.px);
  }

  // Boosts
  __device__ Vec4 Boost(const Vec4& v) const {
    double rsq = M();
    double v0 = (E * v.E - px * v.px - py * v.py - pz * v.pz) / rsq;
    double c1 = (v.E + v0) / (rsq + E);
    return Vec4(v0, v.px - c1 * px, v.py - c1 * py, v.pz - c1 * pz);
  }

  __device__ Vec4 BoostBack(const Vec4& v) const {
    double rsq = M();
    double v0 = (E * v.E + px * v.px + py * v.py + pz * v.pz) / rsq;
    double c1 = (v.E + v0) / (rsq + E);
    return Vec4(v0, v.px + c1 * px, v.py + c1 * py, v.pz + c1 * pz);
  }
};

#endif  // VEC4_CUH_
