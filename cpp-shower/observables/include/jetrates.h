#ifndef JETRATES_H_
#define JETRATES_H_

#include "event.h"

// Jet Rates using the Durham Algorithm
double Yij(const Vec4& p, const Vec4& q);
void Cluster(Event& ev);

#endif  // JETRATES_H_