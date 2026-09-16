#ifndef log_scaled_mom_h_
#define log_scaled_mom_h_

#include "event.h"
#include "histogram.h"

// ---------------------------------------------------------------------------
// Log scaled momentum, xi = -ln(p / p_beam), filled straight into a histogram
// as there is one entry per charged particle

void fill_log_scaled_mom(const event& ev, histo1d& hist, double weight);

#endif  // log_scaled_mom_h_
