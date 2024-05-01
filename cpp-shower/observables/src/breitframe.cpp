#include "breitframe.h"

void CalculateBFObs(Event& ev) {
  // For DIS we know that the electrons are elements 1 and 3 of the event, and
  // the initial state quark is element 2. We can ignore them.

  if (!ev.GetValidity() || ev.GetPartonSize() < 3) {
    return;
  }

  Vec4 moms[maxPartons];
  for (int i = 0; i < ev.GetSize(); ++i) {
    if (i == 0 || i == 1 || i == 3) {
      continue;
    }

    // Only consider partons in the Hemisphere of the Scattered Quark, or
    // the hemisphere along the -z axis.
    if (ev.GetParton(i).GetMom()[3] > 0) {
      continue;
    }

    moms[i - 2] = ev.GetParton(i).GetMom();
  }

  Vec4 z_axis = Vec4(0., 0., 0., 1.);

  double thrust = 0.;
  double e_tot = 0.;
  for (int i = 0; i < ev.GetPartonSize() - 3; ++i) {
    thrust += moms[i].Dot(z_axis);
    e_tot += moms[i].P();  // == E in massless limit
  }

  thrust /= e_tot;

  // ev.SetBF(thrust, jetmas, broadn, cparam);
}