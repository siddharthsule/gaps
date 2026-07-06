#ifndef hadronisation_h_
#define hadronisation_h_

// event and qcd includes all the necessary headers
#include "cluster.h"
#include "event.h"
#include "hadrons.h"
#include "qcd.h"

class hadronisation {
  /**
   * @class hadronisation
   * @brief the hadronisation class
   *
   * Here's to new beginnings! This class will be responsible for a simple
   * cluster hadronisation model
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

  // Kallen Function for Decay Kinematics (isotropic)
  void kallen(vec4 p0, double m1, double m2, vec4& p1, vec4& p2, double rho_1,
              double rho_2) const;

  // Kallen Function for Decay Kinematics (collinear along axis)
  void kallen(vec4 p0, double m1, double m2, vec4& p1, vec4& p2,
              vec4 axis) const;

  // Constituent Reshuffling
  double f_reshuffling(double k, double* masses, double ecms, event& ev) const;
  void constituent_reshuffling(event& ev) const;

  // Flavour Selection for Gluon Splitting, Fission and Decay
  int select_qq_flavour(double rho, bool include_diquarks = false,
                        double rho_diq = 0.0) const;

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
    constituent_reshuffling(ev);
    force_gluons_to_split(ev);
    cluster_list cl;
    form_clusters(ev, cl);
    fission_clusters(ev, cl);
    decay_clusters(ev, cl);
  }
};

#endif  // hadronisation_h_
