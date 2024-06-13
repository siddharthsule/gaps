#ifndef vec4_h_
#define vec4_h_

// base class, with all the important definitions
#include "base.h"

/**
 * four momenta
 * ------------
 *
 * this file contains the definition of the four momenta class, which is used to
 * represent the four-momentum of particles in the event. it is a simple class
 * with the four attributes (e, px, py, pz) and some basic operations that can
 * be performed with them.
 */

class vec4 {
 private:
  double e, px, py, pz;

 public:
  // constructor - define key attributes energy and momentum
  vec4(double e = 0., double px = 0., double py = 0., double pz = 0.)
      : e(e), px(px), py(py), pz(pz) {}

  // get method to obtain attribute value
  double operator[](int i) const {
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
        return 0;
    }
  }

  // print a column vector with the attributes
  friend std::ostream& operator<<(std::ostream& os, const vec4& v) {
    os << "(" << v.e << "," << v.px << "," << v.py << "," << v.pz << ")";
    return os;
  }

  // simple mathematics with four vectors
  vec4 operator+(const vec4& v) const {
    return vec4(e + v.e, px + v.px, py + v.py, pz + v.pz);
  }

  vec4 operator-() const { return vec4(-e, -px, -py, -pz); }

  vec4 operator-(const vec4& v) const {
    return vec4(e - v.e, px - v.px, py - v.py, pz - v.pz);
  }

  // multiplication (and dot product)
  double operator*(const vec4& v) const {
    return e * v.e - px * v.px - py * v.py - pz * v.pz;
  }

  vec4 operator*(double v) const { return vec4(e * v, px * v, py * v, pz * v); }

  // division
  vec4 operator/(double v) const { return vec4(e / v, px / v, py / v, pz / v); }

  // magnitude of the vector
  double m2() const { return (*this) * (*this); }

  double m() const {
    double m2_val = m2();
    return m2_val > 0 ? sqrt(m2_val) : 0;
  }

  // 3 momenta
  double p2() const { return px * px + py * py + pz * pz; }

  double p() const {
    double p2_val = p2();
    return p2_val > 0 ? sqrt(p2_val) : 0;
  }

  // transverse momenta
  double pt2() const { return px * px + py * py; }

  double pt() const {
    double pt2_val = pt2();
    return pt2_val > 0 ? sqrt(pt2_val) : 0;
  }

  // angles
  double theta() const {
    double p_val = p();
    return p_val != 0 ? acos(pz / p_val) : 0;
  }

  double phi() const {
    if (px == 0 && py == 0) {
      return 0.;
    } else {
      return atan2(py, px);
    }
  }

  double rapidity() const {
    double denominator = (e - pz);
    return denominator != 0 ? 0.5 * log((e + pz) / denominator) : 0;
  }

  double eta() const {
    double theta_val = theta();
    return -log(tan(theta_val / 2.));
  }

  // three momenta dot and cross product
  double dot(const vec4& v) const { return px * v.px + py * v.py + pz * v.pz; }

  vec4 cross(const vec4& v) const {
    return vec4(0., py * v.pz - pz * v.py, pz * v.px - px * v.pz,
                px * v.py - py * v.px);
  }

  // boosts
  vec4 dipole_boost(const vec4& v) const {
    double rsq = m();
    double v0 = (e * v.e - px * v.px - py * v.py - pz * v.pz) / rsq;
    double c1 = (v.e + v0) / (rsq + e);
    return vec4(v0, v.px - c1 * px, v.py - c1 * py, v.pz - c1 * pz);
  }

  vec4 dipole_boost_back(const vec4& v) const {
    double rsq = m();
    double v0 = (e * v.e + px * v.px + py * v.py + pz * v.pz) / rsq;
    double c1 = (v.e + v0) / (rsq + e);
    return vec4(v0, v.px + c1 * px, v.py + c1 * py, v.pz + c1 * pz);
  }
};

#endif  // vec4_h_
