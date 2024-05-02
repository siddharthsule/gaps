#ifndef vec4_cuh_
#define vec4_cuh_

// base class, with all the important definitions
#include "base.cuh"

/**
 * four momenta
 * ------------
 *
 * this file contains the definition of the four momenta class, which is used to
 * represent the four-momentum of particles in the event. it is a simple class
 * with the four attributes (e, px, py, pz) and some basic operations that can
 * be performed with them.
 *
 *
 * why is this header only?
 * ------------------------
 *
 * when you declare a __device__ function in a .cuh (header) file and define it
 * in a .cu file, that function's definition is only available to the .cu file
 * in which it is defined. this is because __device__ functions are compiled by
 * nvcc into the device code, and unlike host functions, they do not have
 * external linkage that allows them to be seen or linked across different .cu
 * files after individual compilation.
 */

class vec4 {
 private:
  double e, px, py, pz;

 public:
  // constructor - define key attributes energy and momentum
  // used in me and out [host + device]
  __device__ vec4(double e = 0., double px = 0., double py = 0., double pz = 0.)
      : e(e), px(px), py(py), pz(pz) {}

  // get method to obtain attribute value
  __device__ double operator[](int i) const {
    switch (i) {
      case 0:
        return e;
      case 1:
        return px;
      case 2:
        return py;
      case 3:
        return pz;
      default:
        // cuda does not support exceptions, so we just return 0
        return 0;
    }
  }

  // print a column vector with the attributes
  friend std::ostream& operator<<(std::ostream& os, const vec4& v) {
    os << "(" << v.e << "," << v.px << "," << v.py << "," << v.pz << ")";
    return os;
  }

  // simple mathematics with four vectors
  __device__ vec4 operator+(const vec4& v) const {
    return vec4(e + v.e, px + v.px, py + v.py, pz + v.pz);
  }

  __device__ vec4 operator-() const { return vec4(-e, -px, -py, -pz); }

  __device__ vec4 operator-(const vec4& v) const {
    return vec4(e - v.e, px - v.px, py - v.py, pz - v.pz);
  }

  // multiplication (and dot product)
  __device__ double operator*(const vec4& v) const {
    return e * v.e - px * v.px - py * v.py - pz * v.pz;
  }

  __device__ vec4 operator*(double v) const {
    return vec4(e * v, px * v, py * v, pz * v);
  }

  // division
  __device__ vec4 operator/(double v) const {
    return vec4(e / v, px / v, py / v, pz / v);
  }

  // magnitude of the vector
  __device__ double m2() const { return (*this) * (*this); }

  __device__ double m() const {
    double m2_val = m2();
    return m2_val > 0 ? sqrt(m2_val) : 0;
  }

  // 3 momenta
  __device__ double p2() const { return px * px + py * py + pz * pz; }

  __device__ double p() const {
    double p2_val = p2();
    return p2_val > 0 ? sqrt(p2_val) : 0;
  }

  // transverse momenta
  __device__ double pt2() const { return px * px + py * py; }

  __device__ double pt() const {
    double pt2_val = pt2();
    return pt2_val > 0 ? sqrt(pt2_val) : 0;
  }

  // angles
  __device__ double theta() const {
    double p_val = p();
    return p_val != 0 ? acos(pz / p_val) : 0;
  }

  __device__ double phi() const {
    if (px == 0 && py == 0) {
      return 0.;
    } else {
      return atan2(py, px);
    }
  }

  __device__ double rapidity() const {
    double denominator = (e - pz);
    return denominator != 0 ? 0.5 * log((e + pz) / denominator) : 0;
  }

  __device__ double eta() const {
    double theta_val = theta();
    return -log(tan(theta_val / 2.));
  }

  // three momenta dot and cross product
  __device__ double dot(const vec4& v) const {
    return px * v.px + py * v.py + pz * v.pz;
  }

  __device__ vec4 cross(const vec4& v) const {
    return vec4(0., py * v.pz - pz * v.py, pz * v.px - px * v.pz,
                px * v.py - py * v.px);
  }

  // boosts
  __device__ vec4 boost(const vec4& v) const {
    double rsq = m();
    double v0 = (e * v.e - px * v.px - py * v.py - pz * v.pz) / rsq;
    double c1 = (v.e + v0) / (rsq + e);
    return vec4(v0, v.px - c1 * px, v.py - c1 * py, v.pz - c1 * pz);
  }

  __device__ vec4 boost_back(const vec4& v) const {
    double rsq = m();
    double v0 = (e * v.e + px * v.px + py * v.py + pz * v.pz) / rsq;
    double c1 = (v.e + v0) / (rsq + e);
    return vec4(v0, v.px + c1 * px, v.py + c1 * py, v.pz + c1 * pz);
  }
};

#endif  // vec4_cuh_
