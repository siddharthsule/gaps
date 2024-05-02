#include "dalitz.h"

// dalitz plot

void calculate_dalitz(event& ev) {
  if (!ev.get_validity() || ev.get_parton_size() != 3) {
    return;
  }

  // obtain energy from incoming partons
  double e = abs(ev.get_parton(0).get_mom()[0] + ev.get_parton(1).get_mom()[0]);

  // by default, element 2 is quark and 3 is antiquark
  // i.e. emission will be element 4
  vec4 p1 = ev.get_parton(2).get_mom();
  vec4 p2 = ev.get_parton(3).get_mom();

  // calculate x1 and x2
  double x1 = 2 * p1.p() / e;
  double x2 = 2 * p2.p() / e;

  ev.set_dalitz(x1, x2);
}