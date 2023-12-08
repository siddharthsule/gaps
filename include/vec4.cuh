#ifndef VEC4_CUH_
#define VEC4_CUH_

#include <cuda_runtime.h>

#include <cmath>
#include <iostream>

class Vec4 {
public:

    // Constructor - Define key attributes Energy and Momentum
    // Used in ME and out [HOST + DEVICE]
    __host__ __device__ Vec4(double E=0., double px=0., double py=0., double pz=0.) : E(E), px(px), py(py), pz(pz) {}

    // All others are not used in ME, so HOST only

    // Print a Column Vector with the attributes
    friend std::ostream& operator<<(std::ostream& os, const Vec4& v);

    // Get Method to Obtain Attribute Value
    double operator[](int i) const;

    // Simple Mathematics with Four vectors
    Vec4 operator+(const Vec4& v) const;
    Vec4 operator-(const Vec4& v) const;
    Vec4 operator-() const;

    // Multiplication (and Dot Product)
    double operator*(const Vec4& v) const;
    Vec4 operator*(double v) const;

    Vec4 operator/(double v) const;

    // Magnitude of the Vector
    double M2() const;
    double M() const;

    // Three Momenta
    double P2() const;
    double P() const;

    // Transverse Momenta
    double PT2() const;
    double PT() const;

    // Angles
    double Theta() const;
    double Phi() const;
    double Rapidity() const;
    double Eta() const;

    // Cross Product
    Vec4 Cross(const Vec4& v) const;

    // Boost
    Vec4 Boost(const Vec4& v) const;
    Vec4 BoostBack(const Vec4& v) const;

private:
    double E, px, py, pz;
};

#endif // VEC4_CUH_
