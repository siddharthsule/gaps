#ifndef vec4_cuh_
#define vec4_cuh_

// base class, with all the important definitions
#include "base.cuh"

/**
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
  /**
   * @class vec4
   * @brief The four momentum object.
   *
   * This class is used to store the data for a four momentum. It contains the
   * energy and the three momentum components. It is a simple class with the
   * four attributes (e, px, py, pz) and some basic operations that can be
   * performed with them.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables
  double e, px, py, pz;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  __device__ vec4(double e = 0., double px = 0., double py = 0., double pz = 0.)
      : e(e), px(px), py(py), pz(pz) {}

  // ---------------------------------------------------------------------------
  // getters

  // get method to obtain attribute value
  __host__ __device__ double operator[](int i) const {
    /**
     * @brief get the component of the four vector
     *
     * @param i the index of the component
     * @return the component of the four vector
     */

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
    /**
     * @brief print the four vector
     *
     * @param os the output stream
     * @param v the four vector
     * @return the output stream
     */

    os << "(" << v.e << "," << v.px << "," << v.py << "," << v.pz << ")";
    return os;
  }

  // simple mathematics with four vectors
  __device__ vec4 operator+(const vec4& v) const {
    /**
     * @brief addition of two four vectors
     *
     * @param v the four vector to add
     * @return the sum of the two four vectors
     */

    return vec4(e + v.e, px + v.px, py + v.py, pz + v.pz);
  }

  __device__ vec4 operator-() const {
    /**
     * @brief negation of the four vector
     *
     * @return the negation of the four vector
     */

    return vec4(-e, -px, -py, -pz);
  }

  __device__ vec4 operator-(const vec4& v) const {
    /**
     * @brief subtraction of two four vectors
     *
     * @param v the four vector to subtract
     * @return the difference of the two four vectors
     */

    return vec4(e - v.e, px - v.px, py - v.py, pz - v.pz);
  }

  __device__ double operator*(const vec4& v) const {
    /**
     * @brief dot product of two four vectors
     *
     * @param v the four vector to dot with
     * @return the dot product of the two four vectors
     */

    return e * v.e - px * v.px - py * v.py - pz * v.pz;
  }

  __device__ vec4 operator*(double v) const {
    /**
     * @brief multiplication of a four vector by a scalar
     *
     * @param v the scalar to multiply by
     * @return the four vector multiplied by the scalar
     */

    return vec4(e * v, px * v, py * v, pz * v);
  }

  // division
  __device__ vec4 operator/(double v) const {
    /**
     * @brief division of a four vector by a scalar
     *
     * @param v the scalar to divide by
     * @return the four vector divided by the scalar
     */

    return vec4(e / v, px / v, py / v, pz / v);
  }

  // magnitude of the vector
  __device__ double m2() const {
    /**
     * @brief get the square of the magnitude of the four vector
     *
     * @return the square of the magnitude of the four vector
     */

    return (*this) * (*this);
  }

  __device__ double m() const {
    /**
     * @brief get the magnitude of the four vector
     *
     * @return the magnitude of the four vector
     */

    double m2_val = m2();
    return m2_val > 0 ? sqrt(m2_val) : 0;
  }

  // 3 momenta
  __device__ double p2() const {
    /**
     * @brief get the square of the three momentum
     *
     * @return the square of the three momentum
     */

    return px * px + py * py + pz * pz;
  }

  __device__ double p() const {
    /**
     * @brief get the three momentum
     *
     * @return the three momentum
     */

    double p2_val = p2();
    return p2_val > 0 ? sqrt(p2_val) : 0;
  }

  __device__ double pt2() const {
    /**
     * @brief get the square of the transverse momentum
     *
     * @return the square of the transverse momentum
     */

    return px * px + py * py;
  }

  __device__ double pt() const {
    /**
     * @brief get the transverse momentum
     *
     * @return the transverse momentum
     */

    double pt2_val = pt2();
    return pt2_val > 0 ? sqrt(pt2_val) : 0;
  }

  __device__ double theta() const {
    /**
     * @brief get the polar angle of the three momentum
     *
     * @return the polar angle of the three momentum
     */

    double p_val = p();
    return p_val != 0 ? acos(pz / p_val) : 0;
  }

  __device__ double phi() const {
    /**
     * @brief get the azimuthal angle of the three momentum
     *
     * @return the azimuthal angle of the three momentum
     */

    if (px == 0 && py == 0) {
      return 0.;
    } else {
      return atan2(py, px);
    }
  }

  __device__ double delta_phi(const vec4& v) const {
    /**
     * @brief get the difference in azimuthal angle of the three momentum
     *
     * @param v the four vector to compare with
     * @return the difference in azimuthal angle of the three momentum
     */

    double dphi = phi() - v.phi();
    return dphi > M_PI    ? dphi - 2 * M_PI
           : dphi < -M_PI ? dphi + 2 * M_PI
                          : dphi;
  }

  __device__ double rapidity() const {
    /**
     * @brief get the rapidity of the four vector
     *
     * @return the rapidity of the four vector
     */

    double denominator = (e - pz);
    return denominator != 0 ? 0.5 * log((e + pz) / denominator) : 0;
  }

  __device__ double eta() const {
    /**
     * @brief get the pseudorapidity of the four vector
     *
     * @return the pseudorapidity of the four vector
     */

    double theta_val = theta();
    return -log(tan(theta_val / 2.));
  }

  __device__ double dot(const vec4& v) const {
    /**
     * @brief dot product of the three momenta
     *
     * @param v the three momentum to dot with
     * @return the dot product of the three momenta
     */

    return px * v.px + py * v.py + pz * v.pz;
  }

  __device__ vec4 cross(const vec4& v) const {
    /**
     * @brief cross product of the three momenta
     *
     * @param v the three momentum to cross with
     * @return the cross product of the three momenta
     */

    return vec4(0., py * v.pz - pz * v.py, pz * v.px - px * v.pz,
                px * v.py - py * v.px);
  }

  __device__ vec4 boost(const vec4& v) const {
    /**
     * @brief boost the four vector v to the rest frame of the four vector
     *
     * @param v the four vector to boost
     * @return the boosted four vector
     */

    double rsq = m();
    double v0 = (e * v.e - px * v.px - py * v.py - pz * v.pz) / rsq;
    double c1 = (v.e + v0) / (rsq + e);
    return vec4(v0, v.px - c1 * px, v.py - c1 * py, v.pz - c1 * pz);
  }

  __device__ vec4 boost_back(const vec4& v) const {
    /**
     * @brief boost the four vector v back from the rest frame of the four
     * vector
     *
     * @param v the four vector to boost back
     * @return the boosted four vector
     */

    double rsq = m();
    double v0 = (e * v.e + px * v.px + py * v.py + pz * v.pz) / rsq;
    double c1 = (v.e + v0) / (rsq + e);
    return vec4(v0, v.px + c1 * px, v.py + c1 * py, v.pz + c1 * pz);
  }
};

__device__ inline vec4 operator*(double lhs, const vec4& rhs) {
  return rhs * lhs;
}

#endif  // vec4_cuh_
