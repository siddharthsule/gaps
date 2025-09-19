// -*- C++ -*-
/**
 * GAPS - Event Shapes and Jet Rates at LEP
 * ----------------------------------------
 *
 * In this code, we use the final states of our generated events to calculate
 * some interesting event shapes and jet rates. We use the Durham algorithm to
 * cluster the final state particles into jets and calculate the jet rates. In
 * terms of event shapes, we calculate the thrust, heavy and light jet masses,
 * as well as wide and narrow jet broadenings
 *
 * This code is a simplified version of the ALEPH_2004_I636645 analysis, which
 * includes ALL event shapes and jet rates at LEP I and II.
 */

// Standard
#include <cmath>
#include <iostream>

// Rivet
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/Thrust.hh"

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
class GAPS_LEP : public Analysis {
 private:
  Histo1DPtr _h_y_Durham[5];  // Jets with Durham Algorithm
  Histo1DPtr _h_thrust;       // Thrust
  Histo1DPtr _h_thrust_zoom;  // Jet Masses and Broadenings
  Histo1DPtr _h_heavy_jet_mass;
  Histo1DPtr _h_light_jet_mass;
  Histo1DPtr _h_wide_jet_broadening;
  Histo1DPtr _h_narrow_jet_broadening;
  // Histo2DPtr _h2_dalitz;  // Dalitz Plot for single emission cases

 public:
  /// Constructor
  RIVET_DEFAULT_ANALYSIS_CTOR(GAPS_LEP);

  /// @name Analysis methods
  /// @{

  /// Book histograms and initialise projections before the run
  void init() {
    // Initialise and register projections

    // The basic final-state projection:
    // all final-state particles within
    // the given eta acceptance
    const FinalState fs;
    declare(fs, "FS");

    FastJets durhamjets(fs, FastJets::DURHAM, 0.7, JetAlg::Muons::NONE,
                        JetAlg::Invisibles::NONE);
    declare(durhamjets, "DurhamJets");

    const Thrust thrust(fs);
    declare(thrust, "Thrust");
    declare(Hemispheres(thrust), "Hemispheres");

    // Book histograms
    // specify custom binning
    // book(_h_y_Durham[0], "log10y12", 100, -5., 0.);
    book(_h_y_Durham[1], "log10y23", 100, -4.3, -0.3);
    book(_h_y_Durham[2], "log10y34", 100, -4.3, -0.3);
    book(_h_y_Durham[3], "log10y45", 100, -4.3, -0.3);
    book(_h_y_Durham[4], "log10y56", 100, -4.3, -0.3);
    book(_h_thrust, "tvalue", 100, 0., .5);
    book(_h_thrust_zoom, "tzoomd", 100, 0., .5);
    book(_h_heavy_jet_mass, "hjm", 100, 0., 1.);
    book(_h_light_jet_mass, "ljm", 100, 0., .5);
    book(_h_wide_jet_broadening, "wjb", 100, 0., .5);
    book(_h_narrow_jet_broadening, "njb", 100, 0., .2);
    // book(_h2_dalitz, "dalitz", 100, 0., 1., 100, 0., 1.);
  }

  /// Perform the per-event analysis
  void analyze(const Event& e) {
    // NB: Weight is now Deprecated!

    // Get the Partons if needed for any Analyses
    const FinalState fs = apply<FinalState>(e, "FS");
    Particles partons = fs.particles();

    // Jet Resolutions usng Durham Algorithm
    const FastJets durjet = apply<FastJets>(e, "DurhamJets");
    double log10ynm;
    for (size_t i = 1; i < 5; ++i) {
      double ynm = durjet.clusterSeq()->exclusive_ymerge_max(i + 1);
      if (ynm <= 0.0) {
        _h_y_Durham[i]->fill(-50.0);
      } else {
        log10ynm = log10(ynm);
        _h_y_Durham[i]->fill(log10ynm);
      }
    }

    // Thrust
    const Thrust thrust = apply<Thrust>(e, "Thrust");
    const Vector3 n = thrust.thrustAxis();

    double thr = 1.0 - thrust.thrust();
    if (partons.size() > 2) {
      _h_thrust->fill(thr);
      _h_thrust_zoom->fill(thr);
    } else {
      _h_thrust->fill(-50.);
      _h_thrust_zoom->fill(-50.);
    }

    // Jet Masses and Broadenings
    const Hemispheres& hemi = apply<Hemispheres>(e, "Hemispheres");
    if (partons.size() > 2) {
      _h_heavy_jet_mass->fill(sqrt(hemi.scaledM2high()));
      _h_wide_jet_broadening->fill(hemi.Bmax());

      if (partons.size() > 3) {
        double ljm = sqrt(hemi.scaledM2low());
        double njb = hemi.Bmin();
        if (ljm < 1e-6) ljm = -50.;
        if (njb < 1e-6) njb = -50.;
        _h_light_jet_mass->fill(ljm);
        _h_narrow_jet_broadening->fill(njb);

      } else {
        _h_light_jet_mass->fill(-50.);
        _h_narrow_jet_broadening->fill(-50.);
      }
    } else {
      _h_heavy_jet_mass->fill(-50.);
      _h_wide_jet_broadening->fill(-50.);
      _h_light_jet_mass->fill(-50.);
      _h_narrow_jet_broadening->fill(-50.);
    }

    // // Dalitz Plot
    // if (partons.size() == 3) {
    //   double dal[2] = {0, 0};

    //   Particle q, qbar, g;
    //   bool foundQ(false), foundQBar(false), foundG(false);

    //   for (const Particle& p : partons) {
    //     long id = p.pid();
    //     if ((id != 21) && (id > 0) && !foundQ) {
    //       foundQ = true;
    //       q = p;
    //     } else if ((id != 21) && (id < 0) && !foundQBar) {
    //       foundQBar = true;
    //       qbar = p;
    //     } else if ((id == 21) && !foundG) {
    //       foundG = true;
    //       g = p;
    //     }
    //   }

    //   if (foundQ && foundQBar && foundG) {
    //     // x axis = 2 E_quark / sqrt(s)
    //     // y axis = 2 E_antiquark / sqrt(s)

    //     dal[0] = 2. * q.mom().E() / 91.2;
    //     dal[1] = 2. * qbar.mom().E() / 91.2;

    //     _h2_dalitz->fill(dal[0], dal[1]);
    //   }
    // }
  }

  // Normalise histograms etc., after the run
  void finalize() {
    // n = 0 for y12
    for (size_t n = 1; n < 5; ++n) {
      normalize(_h_y_Durham[n]);
    }

    normalize(_h_thrust);
    normalize(_h_thrust_zoom);
    normalize(_h_heavy_jet_mass);
    normalize(_h_light_jet_mass);
    normalize(_h_wide_jet_broadening);
    normalize(_h_narrow_jet_broadening);
    // normalize(_h2_dalitz);
  }

  /// @}

  /// @name Histograms
  /// @{
  map<string, Histo1DPtr> _h;
  // map<string, Histo2DPtr> _h2;
  map<string, Profile1DPtr> _p;
  map<string, CounterPtr> _c;
  /// @}
};

RIVET_DECLARE_PLUGIN(GAPS_LEP);

}  // namespace Rivet
