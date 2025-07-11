#ifndef eventshapes_h_
#define eventshapes_h_

#include "event.h"

// LEP event shapes
void bubble_sort(vec4* moms, int n);
void calculate_ev_shapes(const event& ev, double* results);

#endif  // eventshapes_h_