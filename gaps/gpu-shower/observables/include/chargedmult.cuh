#ifndef chargedmult_cuh_
#define chargedmult_cuh_

#include "event.cuh"
#include "hadrons.cuh"

// Count charged final-state hadrons for each event and store in results slot.
// Must be called after hadronisation has run; on parton-level events the
// result will be zero (parton PIDs don't appear in the hadron tables).
__global__ void calculate_chargedmult(const event* events, double* results,
                                      int result_slot, int n);

#endif  // chargedmult_cuh_
