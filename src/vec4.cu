#include "vec4.cuh"

// Get Method to Obtain Attribute Value
double Vec4::operator[](int i) const {
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
std::ostream& operator<<(std::ostream& os, const Vec4& v) {
    os << "(" << v.E << "," << v.px << "," << v.py << "," << v.pz << ")";
    return os;
}

// Simple Mathematics with Four vectors
Vec4 Vec4::operator+(const Vec4& v) const {
    return Vec4(E+v.E, px+v.px, py+v.py, pz+v.pz);
}

Vec4 Vec4::operator-() const {
    return Vec4(-E, -px, -py, -pz);
}

Vec4 Vec4::operator-(const Vec4& v) const {
    return Vec4(E-v.E, px-v.px, py-v.py, pz-v.pz);
}

// Multiplication (and Dot Product)
double Vec4::operator*(const Vec4& v) const {
    return E*v.E-px*v.px-py*v.py-pz*v.pz;
}

Vec4 Vec4::operator*(double v) const {
    return Vec4(E*v, px*v, py*v, pz*v);
}

// Division
Vec4 Vec4::operator/(double v) const {
    return Vec4(E/v, px/v, py/v, pz/v);
}

// Magnitude of the Vector
double Vec4::M2() const {
    return (*this)*(*this);
}

double Vec4::M() const {
    double m2 = M2();
    return m2 > 0 ? sqrt(m2) : 0;
}

double Vec4::P2() const {
    return px*px + py*py + pz*pz;
}

double Vec4::P() const {
    double p2 = P2();
    return p2 > 0 ? sqrt(p2) : 0;
}

double Vec4::PT2() const {
    return px*px + py*py;
}

double Vec4::PT() const {
    double pt2 = PT2();
    return pt2 > 0 ? sqrt(pt2) : 0;
}

double Vec4::Theta() const {
    double p = P();
    return p != 0 ? acos(pz/p) : 0;
}

double Vec4::Phi() const {
    if (px == 0 && py == 0) {
        return 0.0;
    } else {
        return atan2(py, px);
    }
}

double Vec4::Rapidity() const {
    double denominator = (E - pz);
    return denominator != 0 ? 0.5 * log((E + pz)/denominator) : 0;
}

double Vec4::Eta() const {
    double theta = Theta();
    return - log(tan(theta/2.));
}

Vec4 Vec4::Cross(const Vec4& v) const {
    return Vec4(0.0,
                py*v.pz - pz*v.py,
                pz*v.px - px*v.pz,
                px*v.py - py*v.px);
}

Vec4 Vec4::Boost(const Vec4& v) const {
    double rsq = M();
    double v0 = (E*v.E - px*v.px - py*v.py - pz*v.pz)/rsq;
    double c1 = (v.E + v0)/(rsq + E);
    return Vec4(v0,
                v.px - c1*px,
                v.py - c1*py,
                v.pz - c1*pz);
}

Vec4 Vec4::BoostBack(const Vec4& v) const {
    double rsq = M();
    double v0 = (E*v.E + px*v.px + py*v.py + pz*v.pz)/rsq;
    double c1 = (v.E + v0)/(rsq + E);
    return Vec4(v0,
                v.px + c1*px,
                v.py + c1*py,
                v.pz + c1*pz);
}
