#include "breitframe.h"

void calculate_bf_obs(event& ev) {
  // for DIS we know that the electrons are elements 1 and 3 of the event, and
  // the initial state quark is element 2. we can ignore them.

  if (!ev.get_validity() || ev.get_parton_size() < 3) {
    return;
  }

  // get the useful final state partons
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

  // define the "z axis" as the axis facing the -z direction
  vec4 z_axis = vec4(0., 0., 0., -1.);

  double e_tot = 0.;
  double thrust = 0.;
  double jetmas = 0.;
  double broadn = 0.;

  vec4 p_with = vec4();

  for (int i = 0; i < ev.get_parton_size() - 3; ++i) {
    // Energy
    e_tot += moms[i].p();  // == e in massless limit

    // Thrust
    thrust += moms[i].dot(z_axis);

    // Jet Mass
    p_with = p_with + moms[i];

    // Broadening
    double mo_para = moms[i].dot(z_axis);
    double mo_perp = (moms[i] - (z_axis * mo_para)).p();
    broadn += mo_perp;
  }

  thrust /= e_tot;
  thrust = 1. - thrust;

  jetmas = fabs(p_with.m2() / (4 * e_tot * e_tot));
  jetmas = sqrt(jetmas);

  broadn /= (2 * e_tot);

  // ev.set_bf(thrust, jetmas, broadn);
}