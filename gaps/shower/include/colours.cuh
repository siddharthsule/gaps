#ifndef SHOWER_COLOURS_CUH_
#define SHOWER_COLOURS_CUH_

#include "qcd.cuh"

// Colours
__device__ void MakeColours(Event &ev, int *coli, int *colj, const int flavs[3],
                            const int colij[2], const int colk[2],
                            const double rand);

#endif // SHOWER_COLOURS_CUH_