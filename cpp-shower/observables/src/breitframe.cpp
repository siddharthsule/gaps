#include "breitframe.h"

void calculate_bf_obs(event& ev) {
  // for dis we know that the electrons are elements 1 and 3 of the event, and
  // the initial state quark is element 2. we can ignore them.

  if (!ev.get_validity() || ev.get_parton_size() < 3) {
    return;
  }

  vec4 moms[max_partons];
  for (int i = 0; i < ev.get_size(); ++i) {
    if (i == 0 || i == 1 || i == 3) {
      continue;
    }

    // only consider partons in the hemisphere of the scattered quark, or
    // the hemisphere along the -z axis.
    if (ev.get_parton(i).get_mom()[3] > 0) {
      continue;
    }

    moms[i - 2] = ev.get_parton(i).get_mom();
  }

  vec4 z_axis = vec4(0., 0., 0., 1.);

  double thrust = 0.;
  double e_tot = 0.;
  for (int i = 0; i < ev.get_parton_size() - 3; ++i) {
    thrust += moms[i].dot(z_axis);
    e_tot += moms[i].p();  // == e in massless limit
  }

  thrust /= e_tot;

  // ev.set_bf(thrust, jetmas, broadn, cparam);
}