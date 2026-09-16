#ifndef hadronisation_h_
#define hadronisation_h_

// event and qcd includes all the necessary headers
#include "cluster.h"
#include "decay_kinematics.h"
#include "event.h"
#include "qcd.h"
#include "utilities.h"

class hadronisation {
  /**
   * @class hadronisation
   * @brief the cluster hadronisation algorithm
   */

 private:
  // -------------------------------------------------------------------------
  // member variables

  // Constituent masses [d, u, s, c, b] in GeV
  double const_mass[5] = {0.33, 0.33, 0.45, mc, mb};  // GeV
  double gluon_mass = 0.95;                           // GeV

  // Cluster fission parameters (light, charm, bottom)
  double clmax[3] = {3.53, 3.95, 3.76};   // max cluster mass (GeV)
  double clpow[3] = {1.85, 2.56, 0.55};   // power in threshold condition
  double psplit[3] = {0.91, 0.99, 0.63};  // power in mass sampling distribution

  // Probability weights for flavours (d, u, s, diquark)
  double pwt[4] = {1.0, 1.0, 0.37, 0.33};

 public:
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
  // member functions

  // Constituent Reshuffling
  double f_reshuffling(double k, double* masses, double ecms, event& ev) const;
  void constituent_reshuffling(event& ev) const;

  // Gluon Splitting
  void force_gluons_to_split(event& ev) const;

  // Cluster Formation
  void form_clusters(event& ev, cluster_list& cl) const;

  // Cluster Fission
  void fission_clusters(event& ev, cluster_list& cl) const;

  // Cluster Decay
  void decay_clusters(event& ev, cluster_list& cl) const;

  // Wrapper to run the full hadronisation sequence
  void run(event& ev) const {
    /**
     * @brief hadronise one event, from partons to hadrons.
     *
     * Every step skips an overflowed event; validate() drops it later.
     */

    // Constituent Mass Reshuffling
    constituent_reshuffling(ev);

    // Forced Gluon Splitting
    force_gluons_to_split(ev);

    // Cluster Formation
    cluster_list cl;
    form_clusters(ev, cl);

    // Cluster Fission
    fission_clusters(ev, cl);

    // Cluster Decay
    decay_clusters(ev, cl);
  }
};

#endif  // hadronisation_h_
