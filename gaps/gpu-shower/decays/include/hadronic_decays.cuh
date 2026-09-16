#ifndef hadronic_decays_cuh_
#define hadronic_decays_cuh_

#include "database.cuh"
#include "decay_kinematics.cuh"
#include "event.cuh"
#include "interface.cuh"

class hadronic_decays {
  /**
   * @class hadronic_decays
   * @brief Decay the hadrons cluster decay leaves until the event is final.
   *
   * States with ctau >= ctau_max are final; at 10 mm that leaves pi+-, K+-,
   * K_S, K_L, p, n and the weakly decaying light baryons.
   */

 public:
  // ---------------------------------------------------------------------------
  // member variables
  // Public so that __global__ kernels (which are not class members) can read
  // them via the dec* pointer passed at launch time.

  // Proper decay length at or above which a state is a final state particle
  double ctau_max = 10.;  // mm

  // Cap on the passes over the event, so a malformed table cannot loop forever
  int max_generations = 100;

  // ---------------------------------------------------------------------------
  // constructor

  hadronic_decays() = default;

  explicit hadronic_decays(double ctau_max_in) : ctau_max(ctau_max_in) {}

  // ---------------------------------------------------------------------------
  // member functions

  // Is this a final state particle?
  __device__ bool is_stable(int pid) const;

  // Decay one state, returning false where it was left alone
  __device__ bool decay_one(event& ev, int i) const;

  // Decay every unstable particle, and every unstable child of those decays
  __device__ void run(event& ev) const;
};

// -----------------------------------------------------------------------------
// Kernel — one thread per event

__global__ void hadronic_decays_kernel(event* events, hadronic_decays* dec,
                                       int n);

// -----------------------------------------------------------------------------
// Host wrapper — allocates the decay object on device and launches the kernel

void run_decays(thrust::device_vector<event>& dv_events, const params& p,
                int blocks);

#endif  // hadronic_decays_cuh_
