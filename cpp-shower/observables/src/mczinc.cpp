#include "mczinc.h"

void calculate_zboson_obs(event& ev) {
  // for LHC we know that the electrons are elements 2 and 3 of the event. We
  // can create the z boson momentum using the invariant mass of the two (?)

  if (!ev.get_validity() || ev.get_parton_size() < 3) {
    return;
  }

  // calculate the z boson momentum
  vec4 z_mom = ev.get_parton(2).get_mom() + ev.get_parton(3).get_mom();

  // ev.set_z_mass(z_mom.m());
  // ev.set_z_pt(z_mom.pt());
  // ev.set_z_phi(z_mom.phi());
  // ev.set_z_eta(z_mom.eta());
}

void calculate_lepton_obs(event& ev) {
  // for LHC we know that the electrons are elements 2 and 3 of the event. We
  // just need to get their momenta and measure the pt and eta.

  if (!ev.get_validity() || ev.get_parton_size() < 3) {
    return;
  }

  double lep_pt = 0.;
  double lep_eta = 0.;

  for (int i = 2; i < 4; i++) {
    lep_pt += ev.get_parton(i).get_mom().pt();
    lep_eta += ev.get_parton(i).get_mom().eta();
  }

  // ev.set_lep_pt(lep_pt);
  // ev.set_lep_eta(lep_eta);
}