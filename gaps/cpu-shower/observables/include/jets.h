#ifndef jets_h_
#define jets_h_

#include "event.h"

// jet rates using the durham algorithm
double yij(const vec4& p, const vec4& q);
void cluster_durham(const event& ev, double* results);

// jet clustering using the generalized kt algorithm
void bubble_sort_pt(vec4* moms, int n);
double dR2(const vec4& p, const vec4& q);
double dij(const vec4& p, const vec4& q);
void cluster_genkt(const event& ev, double* results);

#endif  // jets_h_