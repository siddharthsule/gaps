#ifndef jets_h_
#define jets_h_

#include "event.h"

// jet rates using the durham algorithm
double yij(const vec4& p, const vec4& q);
void cluster_durham(const event& ev, double* results);

#endif  // jets_h_