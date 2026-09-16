#ifndef hadronic_decays_h_
#define hadronic_decays_h_

#include "database.h"
#include "decay_kinematics.h"
#include "event.h"

class hadronic_decays {
  /**
   * @class hadronic_decays
   * @brief Decay the hadrons cluster decay leaves until the event is final.
   *
   * States with ctau >= ctau_max are final; at 10 mm that leaves pi+-, K+-,
   * K_S, K_L, p, n and the weakly decaying light baryons.
   */

 private:
  // ---------------------------------------------------------------------------
  // member variables

  // Proper decay length at or above which a state is a final state particle
  double ctau_max = 10.;  // mm

  // Cap on the passes over the event, so a malformed table cannot loop forever
  int max_generations = 100;

 public:
  // ---------------------------------------------------------------------------
  // constructor

  hadronic_decays() = default;

  explicit hadronic_decays(double ctau_max_in) : ctau_max(ctau_max_in) {}

  // ---------------------------------------------------------------------------
  // member functions

  // Is this a final state particle?
  bool is_stable(int pid) const;

  // Decay one state, returning false where it was left alone
  bool decay_one(event& ev, int i) const;

  // Decay every unstable particle, and every unstable child of those decays
  void run(event& ev) const;
};

#endif  // hadronic_decays_h_
