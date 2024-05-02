#ifndef jetrates_h_
#define jetrates_h_

#include "event.h"

// jet rates using the durham algorithm
double yij(const vec4& p, const vec4& q);
void cluster(event& ev);

#endif  // jetrates_h_