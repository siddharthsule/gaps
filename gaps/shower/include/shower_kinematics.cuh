#ifndef SHOWER_KINEMATICS_CUH_
#define SHOWER_KINEMATICS_CUH_

#include "qcd.cuh"

// Kinematics
__device__ void MakeKinematics(Vec4 *kinematics, const double z, const double y,
                               const double phi, const Vec4 pijt,
                               const Vec4 pkt);

#endif // SHOWER_KINEMATICS_CUH_