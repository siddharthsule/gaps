// -*- C++ -*-
/**
 * GAPS - Z Boson Observables and Jet Rates at LHC
 * -----------------------------------------------
 *
 * In this code, we focus on the Z boson production and its decay products.
 * we also study the behaviour of the emitted partons. We explicitly remove
 * the beam remnants before doing the analysis too. The observables we look
 * at are the Z pT, phi, y, mass, lepton and jet pT and eta, as well as the
 * angular correlations between the Z boson and the jets.
 *
 * This code is a simplified version of the MC_ZINC and MC_ZJETS analyses.
 */

// Standard
#include <cmath>
#include <iostream>

// Rivet
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

/**
 * STEPS TO BUILD ANALYSIS:
 * ------------------------
 * - Create a new directory for the analysis
 * - rivet-mkanalysis MyAnalysis
 * - Add your code to MyAnalysis.cc
 * - rivet-build RivetMyAnalysis.so MyAnalysis.cc
 * - export RIVET_ANALYSIS_PATH=$PWD
 * - To Check: rivet --show-analyis MyAnalysis
 */

namespace Rivet {

/// @brief Add a short analysis description here
class GAPS_LHC : public Analysis {
 private:
  // MC_ZINC
  Histo1DPtr _h_Z_mass;
  Histo1DPtr _h_Z_pT;
  Histo1DPtr _h_Z_pT_full;
  Histo1DPtr _h_Z_y;
  Histo1DPtr _h_Z_phi;
  Histo1DPtr _h_lepton_pT;
  Histo1DPtr _h_lepton_eta;

  // MC_ZJETS
  Histo1DPtr _h_jet1_pt;
  Histo1DPtr _h_jet2_pt;
  Histo1DPtr _h_jet3_pt;
  Histo1DPtr _h_jet1_eta;
  Histo1DPtr _h_jet2_eta;
  Histo1DPtr _h_jet3_eta;
  Histo1DPtr _h_num_jets;

  // MIX
  Histo1DPtr _h_deta_z_j1;
  Histo1DPtr _h_dR_z_j1;
  Histo1DPtr _h_dR_j1_j2;
  Histo1DPtr _h_dR_j1_j3;
  Histo1DPtr _h_dR_j2_j3;

  double _z_lo = 60.;
  double _z_hi = 120.;
  double _dR;
  PdgId _lepton;

 public:
  /// Constructor
  RIVET_DEFAULT_ANALYSIS_CTOR(GAPS_LHC);

  /// @name Analysis methods
  /// @{

  /// Book histograms
  void init() {
    // z
    book(_h_Z_mass, "zmass", 100, _z_lo, _z_hi);
    book(_h_Z_pT, "zpt", 100, 0., 100.);
    book(_h_Z_pT_full, "zptfull", logspace(100, 0.5, 1000.));
    book(_h_Z_phi, "zphi", 100, -3.2, 3.2);
    book(_h_Z_y, "zrap", 100, -10., 10.);
    book(_h_lepton_pT, "leptpt", 100, 0., 100.);
    book(_h_lepton_eta, "lepteta", 100, -10., 10.);

    // jets
    book(_h_jet1_pt, "j1pt", 100, 0., 100.);
    book(_h_jet2_pt, "j2pt", 100, 0., 100.);
    book(_h_jet3_pt, "j3pt", 100, 0., 100.);
    book(_h_jet1_eta, "j1eta", 100, -10., 10);
    book(_h_jet2_eta, "j2eta", 100, -10., 10);
    book(_h_jet3_eta, "j3eta", 100, -10., 10);
    book(_h_num_jets, "njets", 10, -0.5, 9.5);

    // Mixing z and zjets
    book(_h_deta_z_j1, "detazj1", 100, -10., 10.);
    book(_h_dR_z_j1, "dRzj1", 100, 0., 10.);
    book(_h_dR_j1_j2, "dRj1j2", 100, 0., 10.);
    book(_h_dR_j1_j3, "dRj1j3", 100, 0., 10.);
    book(_h_dR_j2_j3, "dRj2j3", 100, 0., 10.);

    // Assign a few final states for differnt analyses
    // for Z Boson (-> e+ e-)
    IdentifiedFinalState lepton_fs;
    lepton_fs.acceptIdPair(PID::ELECTRON);
    lepton_fs.acceptId(PID::Z0);
    declare(lepton_fs, "LFS");

    VetoedFinalState parton_fs;
    parton_fs.addVeto(PID::ELECTRON);
    parton_fs.addVeto(PID::POSITRON);
    parton_fs.addVeto(PID::Z0);
    parton_fs.addDecayProductsVeto(82);
    declare(parton_fs, "PFS");

    // ZFinder
    _dR = 0.0;  // Placeholder dR
    _lepton = PID::ELECTRON;
    Cut cut = Cuts::abspid == _lepton;  // Placeholder cut
    ZFinder zfinder(lepton_fs, cut, _lepton, _z_lo * GeV, _z_hi * GeV, _dR,
                    ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
    declare(zfinder, "ZFinder");

    // FastJets
    FastJets jets(parton_fs, FastJets::ANTIKT, 0.4);
    declare(jets, "Jets");
  }

  /// Do the analysis
  void analyze(const Event& e) {
    /**
     * This can also be done using the DiLeptonFinder or ZFinder Class, but
     * putting this version to show how you can code things on your own.
     */

    // Get the final state particles
    const FinalState lepton_fs = apply<FinalState>(e, "LFS");
    Particles leptons = lepton_fs.particles();
    const FinalState parton_fs = apply<FinalState>(e, "PFS");
    Particles partons = parton_fs.particles();

    /**
     * ZFinder Version:
     */
    FourMomentum zmom;
    const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
    if (zfinder.bosons().size() == 1) {
      zmom = zfinder.bosons()[0].momentum();
    } else {
      // Manually check for an on-shell Z boson
      bool found = false;
      for (const Particle& p : leptons) {
        if (p.pid() == PID::Z0) {
          found = true;
          zmom = p.momentum();
          break;
        }
      }
      if (!found) {
        vetoEvent;
      }
    }
    _h_Z_mass->fill(zmom.mass() / GeV);
    _h_Z_pT->fill(zmom.pT() / GeV);
    _h_Z_pT_full->fill(zmom.pT() / GeV);
    _h_Z_y->fill(zmom.rapidity());

    if (zmom.phi() < 1e-6) {
      _h_Z_phi->fill(zmom.phi());
    } else {
      _h_Z_phi->fill(zmom.phi() - PI);
    }

    for (const Particle& l : zfinder.constituents()) {
      _h_lepton_pT->fill(l.pT() / GeV, 0.5);
      _h_lepton_eta->fill(l.eta(), 0.5);
    }

    Cut cut = Cuts::pt > 0. * GeV;

    bool do_clustering = true;

    // Do Jet Clustering
    if (do_clustering) {
      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(cut);

      // Extract jet kinematics for applying cuts
      double jet1_pt = jets.size() > 0 ? jets[0].pT() / GeV : -50.;
      double jet1_eta = jets.size() > 0 ? jets[0].eta() : -50.;
      double jet2_pt = jets.size() > 1 ? jets[1].pT() / GeV : -50.;
      double jet2_eta = jets.size() > 1 ? jets[1].eta() : -50.;
      double jet3_pt = jets.size() > 2 ? jets[2].pT() / GeV : -50.;
      double jet3_eta = jets.size() > 2 ? jets[2].eta() : -50.;

      // Count jets satisfying |eta| < 5 and pt > 5
      int n_good_jets = 0;
      for (size_t i = 0; i < jets.size(); ++i) {
        if (fabs(jets[i].eta()) < 5.0 && jets[i].pT() / GeV > 5.0) {
          n_good_jets++;
        }
      }

      // Jet1 observables with kinematic cuts
      if (jets.size() > 0) {
        // jetpt1 only when |jeteta1|<5
        _h_jet1_pt->fill((fabs(jet1_eta) < 5.0) ? jet1_pt : -50.);
        // jeteta1 only when jetpt1>5
        _h_jet1_eta->fill((jet1_pt > 5.0) ? jet1_eta : -50.);

        // Z-jet1 observables only when |jeteta1|<5 and jetpt1>5
        if (fabs(jet1_eta) < 5.0 && jet1_pt > 5.0) {
          _h_deta_z_j1->fill(zmom.eta() - jets[0].eta());
          _h_dR_z_j1->fill(deltaR(zmom, jets[0].momentum()));
        } else {
          _h_deta_z_j1->fill(-50.);
          _h_dR_z_j1->fill(-50.);
        }
      } else {
        _h_deta_z_j1->fill(-50.);
        _h_dR_z_j1->fill(-50.);
        _h_jet1_pt->fill(-50.);
        _h_jet1_eta->fill(-50.);
      }

      // Jet2 observables with kinematic cuts
      if (jets.size() > 1) {
        // jetpt2 only when |jeteta2|<5
        _h_jet2_pt->fill((fabs(jet2_eta) < 5.0) ? jet2_pt : -50.);
        // jeteta2 only when jetpt2>5
        _h_jet2_eta->fill((jet2_pt > 5.0) ? jet2_eta : -50.);

        // dR(j1,j2) only when both jets pass cuts
        if (fabs(jet1_eta) < 5.0 && jet1_pt > 5.0 && fabs(jet2_eta) < 5.0 &&
            jet2_pt > 5.0) {
          _h_dR_j1_j2->fill(deltaR(jets[0].momentum(), jets[1].momentum()));
        } else {
          _h_dR_j1_j2->fill(-50.);
        }
      } else {
        _h_jet2_pt->fill(-50.);
        _h_jet2_eta->fill(-50.);
        _h_dR_j1_j2->fill(-50.);
      }

      // Jet3 observables with kinematic cuts
      if (jets.size() > 2) {
        // jetpt3 only when |jeteta3|<5
        _h_jet3_pt->fill((fabs(jet3_eta) < 5.0) ? jet3_pt : -50.);
        // jeteta3 only when jetpt3>5
        _h_jet3_eta->fill((jet3_pt > 5.0) ? jet3_eta : -50.);

        // dR(j1,j3) only when both jets pass cuts
        if (fabs(jet1_eta) < 5.0 && jet1_pt > 5.0 && fabs(jet3_eta) < 5.0 &&
            jet3_pt > 5.0) {
          _h_dR_j1_j3->fill(deltaR(jets[0].momentum(), jets[2].momentum()));
        } else {
          _h_dR_j1_j3->fill(-50.);
        }

        // dR(j2,j3) only when both jets pass cuts
        if (fabs(jet2_eta) < 5.0 && jet2_pt > 5.0 && fabs(jet3_eta) < 5.0 &&
            jet3_pt > 5.0) {
          _h_dR_j2_j3->fill(deltaR(jets[1].momentum(), jets[2].momentum()));
        } else {
          _h_dR_j2_j3->fill(-50.);
        }
      } else {
        _h_jet3_pt->fill(-50.);
        _h_jet3_eta->fill(-50.);
        _h_dR_j1_j3->fill(-50.);
        _h_dR_j2_j3->fill(-50.);
      }

      // Fill good jets count instead of total jet count
      _h_num_jets->fill(n_good_jets);
    }

    // Do Parton Momenta Only
    else {
      const Particles& partonspt = parton_fs.particlesByPt(cut);

      if (partonspt.size() > 0) {
        _h_deta_z_j1->fill(zmom.eta() - partonspt[0].eta());
        _h_dR_z_j1->fill(deltaR(zmom, partonspt[0].momentum()));
        _h_jet1_pt->fill(partonspt[0].pT() / GeV);
        _h_jet1_eta->fill(partonspt[0].eta());
      } else {
        _h_deta_z_j1->fill(-50.);
        _h_dR_z_j1->fill(-50.);
        _h_jet1_pt->fill(-50.);
        _h_jet1_eta->fill(-50.);
      }

      if (partonspt.size() > 1) {
        _h_jet2_pt->fill(partonspt[1].pT() / GeV);
        _h_jet2_eta->fill(partonspt[1].eta());
        _h_dR_j1_j2->fill(
            deltaR(partonspt[0].momentum(), partonspt[1].momentum()));
      } else {
        _h_jet2_pt->fill(-50.);
        _h_jet2_eta->fill(-50.);
        _h_dR_j1_j2->fill(-50.);
      }

      if (partonspt.size() > 2) {
        _h_jet3_pt->fill(partonspt[2].pT() / GeV);
        _h_jet3_eta->fill(partonspt[2].eta());
        _h_dR_j1_j3->fill(
            deltaR(partonspt[0].momentum(), partonspt[2].momentum()));
        _h_dR_j2_j3->fill(
            deltaR(partonspt[1].momentum(), partonspt[2].momentum()));
      } else {
        _h_jet3_pt->fill(-50.);
        _h_jet3_eta->fill(-50.);
        _h_dR_j1_j3->fill(-50.);
        _h_dR_j2_j3->fill(-50.);
      }

      // Exclusive
      _h_num_jets->fill(partonspt.size());
    }
  }

  /// Finalize
  void finalize() {
    // Normalisation
    // const double norm = (crossSection() / picobarn) / sumOfWeights();
    const double norm = 1. / sumOfWeights();

    scale(_h_Z_mass, norm);
    scale(_h_Z_pT, norm);
    scale(_h_Z_pT_full, norm);
    scale(_h_Z_phi, norm);
    scale(_h_Z_y, norm);
    scale(_h_lepton_pT, norm);
    scale(_h_lepton_eta, norm);

    scale(_h_jet1_pt, norm);
    scale(_h_jet2_pt, norm);
    scale(_h_jet3_pt, norm);
    scale(_h_jet1_eta, norm);
    scale(_h_jet2_eta, norm);
    scale(_h_jet3_eta, norm);
    scale(_h_num_jets, norm);

    scale(_h_deta_z_j1, norm);
    scale(_h_dR_z_j1, norm);
    scale(_h_dR_j1_j2, norm);
    scale(_h_dR_j1_j3, norm);
    scale(_h_dR_j2_j3, norm);
  }

  /// @name Histograms
  /// @{
  map<string, Histo1DPtr> _h;
  map<string, Profile1DPtr> _p;
  map<string, CounterPtr> _c;
  /// @}
};

RIVET_DECLARE_PLUGIN(GAPS_LHC);

}  // namespace Rivet
