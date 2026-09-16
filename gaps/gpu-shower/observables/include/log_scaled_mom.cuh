#ifndef log_scaled_mom_cuh_
#define log_scaled_mom_cuh_

#include "event.cuh"
#include "histogram.cuh"

// ---------------------------------------------------------------------------
// Log scaled momentum, xi = -ln(p / p_beam), filled straight into a histogram
// as there is one entry per charged particle

__device__ void fill_log_scaled_mom(const event& ev, histo1d& hist,
                                    double weight);

#endif  // log_scaled_mom_cuh_
