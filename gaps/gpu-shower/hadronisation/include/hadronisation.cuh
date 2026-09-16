#ifndef hadronisation_cuh_
#define hadronisation_cuh_

// event and qcd includes all the necessary headers
#include "cluster.cuh"
#include "decay_kinematics.cuh"
#include "event.cuh"
#include "interface.cuh"
#include "qcd.cuh"
#include "utilities.cuh"

class hadronisation {
  /**
   * @class hadronisation
   * @brief the cluster hadronisation algorithm
   */

 public:
  // -------------------------------------------------------------------------
  // member variables
  // Public so that __global__ kernels (which are not class members) can read
  // them via the had* pointer passed at launch time.

  // Constituent masses [d, u, s, c, b] in GeV
  double const_mass[5] = {0.33, 0.33, 0.45, mc, mb};  // GeV
  double gluon_mass = 0.95;                           // GeV

  // Cluster fission parameters (light, charm, bottom)
  double clmax[3] = {3.53, 3.95, 3.76};   // max cluster mass (GeV)
  double clpow[3] = {1.85, 2.56, 0.55};   // power in threshold condition
  double psplit[3] = {0.91, 0.99, 0.63};  // power in mass sampling distribution

  // Probability weights for flavours (d, u, s, diquark)
  double pwt[4] = {1.0, 1.0, 0.37, 0.33};

  // -------------------------------------------------------------------------
  // constructor

  hadronisation() = default;

  /**
   * @brief construct a hadronisation model with custom fission and flavour
   * parameters, overriding the defaults.
   */
  hadronisation(const double clmax_in[3], const double clpow_in[3],
                const double psplit_in[3], const double pwt_in[4]) {
    // Set the fission parameters
    for (int i = 0; i < 3; i++) {
      clmax[i] = clmax_in[i];
      clpow[i] = clpow_in[i];
      psplit[i] = psplit_in[i];
    }
    // Set the flavour weights
    for (int i = 0; i < 4; i++) {
      pwt[i] = pwt_in[i];
    }
  }

  // -------------------------------------------------------------------------
  // Device helper functions (called from kernels)

  // Constituent Reshuffling
  __device__ double f_reshuffling(double k, double* masses, double ecms,
                                  event& ev) const;

};

// -----------------------------------------------------------------------------
// The hadrons every flavour pair forms, built once on the device before any
// step that reads them. The table itself is file local to decay_clusters.cu;
// only the launch has to be visible here.

__global__ void build_candidate_table();

// -----------------------------------------------------------------------------
// Kernels — one per hadronisation step, one thread per event

// Constituent Reshuffling
__global__ void constituent_reshuffling(event* events, hadronisation* had,
                                        int n);

// Gluon Splitting
__global__ void force_gluons_to_split(event* events, hadronisation* had, int n);

// Cluster Formation
__global__ void form_clusters(event* events, cluster_list* cls,
                              hadronisation* had, int n);

// Cluster Fission
__global__ void fission_clusters(event* events, cluster_list* cls,
                                 hadronisation* had, int n);

// Cluster Decay
__global__ void decay_clusters(event* events, cluster_list* cls,
                               hadronisation* had, int n);

// -----------------------------------------------------------------------------
// Host wrapper — allocates the hadronisation object on device and launches the
// steps in sequence

void run_hadronisation(thrust::device_vector<event>& dv_events,
                       thrust::device_vector<cluster_list>& dv_cls,
                       const params& p, int blocks);

#endif  // hadronisation_cuh_
