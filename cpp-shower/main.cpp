// to measure wall clock time and write to file
#include <chrono>
#include <fstream>

// base components
#include "base.h"

// me
#include "matrix.h"

// shower
#include "shower.h"

// analysis
#include "observables.h"

/**
 * GAPS: C++ Shower for Comparison
 * ------------------------------------
 *
 * This program is a translation of S. Höche's "Introduction to Parton Showers"
 * Python tutorial[1], with added functionality for parallelisation, a Event
 * class and event shape analyses.
 *
 * The purpose of this program is to compare the performance of the C++ and
 * CUDA versions of the shower, and to compare the performance of the C++ with
 * parallelisation and CUDA.
 *
 * [1] https://arxiv.org/abs/1411.4085 and MCNET-CTEQ 2021 Tutorial
 */

// -----------------------------------------------------------------------------

void run_generator(const int& n, const double& e, const std::string& filename) {
  // ---------------------------------------------------------------------------
  // give some information about the simulation

  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "|        GAPS: C++ Shower for Comparison        |" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Process: e+ e- --> q qbar" << std::endl;
  std::cout << "Number of Events: " << n << std::endl;
  std::cout << "Centre of Mass Energy: " << e << " GeV" << std::endl;
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // inititalisation

  std::cout << "Initialising..." << std::endl;
  std::vector<event> events(n);

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating Matrix Elements (C++)..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  matrix xs(asmz, e);

  for (int i = 0; i < n; i++) {
    xs.generate_lo_point(events[i]);  // random seed
    // (same seed option is off in matrix.cpp!)
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // showering

  std::cout << "Showering Partons (C++)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  shower sh;

  for (int i = 0; i < n; i++) {
    sh.run(events[i]);  // random seed
    // (same seed option is off in shower.cpp!)
  }

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_sh = end - start;

  // ---------------------------------------------------------------------------
  // analysis

  std::cout << "Analysing Events (C++)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // remove existing file
  std::remove(filename.c_str());

  analysis an;

  // analyze events (including validation of colour and momentum conservation)
  for (int i = 0; i < n; i++) {
    an.analyze(events[i]);
  }

  // storage
  an.finalize(filename);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_an = end - start;

  // ---------------------------------------------------------------------------
  // results

  double diff = diff_me.count() + diff_sh.count() + diff_an.count();

  std::cout << "" << std::endl;
  std::cout << "EVENT GENERATION COMPLETE" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "ME Time: " << diff_me.count() << " s" << std::endl;
  std::cout << "Sh Time: " << diff_sh.count() << " s" << std::endl;
  std::cout << "An Time: " << diff_an.count() << " s" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Total time: " << diff << " s" << std::endl;
  std::cout << "" << std::endl;

  // open the file in append mode. this will create the file if it doesn't
  // exist.
  std::ofstream outfile("cpp-time.dat", std::ios_base::app);

  // write diff_sh.count() to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // close the file.
  outfile.close();

  std::cout << "Histograms written to " << filename << std::endl;
  std::cout << "Timing data written to cpp-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  // import settings
  int n = argc > 1 ? atoi(argv[1]) : 10000;
  double e = argc > 2 ? atof(argv[2]) : 91.2;

  // if more than max_events, run in batches
  if (n > max_events) {
    std::cout << "-------------------------------------------------"
              << std::endl;
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    std::cout << "Running in batches..." << std::endl;
    std::cout << "Please use rivet-merge to combine runs" << std::endl;

    // split into batches
    int n_batches = n / max_events;
    int n_remainder = n % max_events;
    std::cout << "Number of Batches: " << n_batches << std::endl;
    std::cout << "Size of Remainder: " << n_remainder << std::endl;

    // run in batches
    for (int i = 0; i < n_batches; i++) {
      std::string filename = "cpp-" + std::to_string(i) + ".yoda";
      run_generator(max_events, e, filename);
    }

    // run remainder
    if (n_remainder > 0) {
      std::string filename = "cpp-" + std::to_string(n_batches) + ".yoda";
      run_generator(n_remainder, e, filename);
    }
  } else {
    run_generator(n, e, "cpp.yoda");
  }

  return 0;
}
// -----------------------------------------------------------------------------
