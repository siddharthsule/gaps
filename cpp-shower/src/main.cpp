// To Measure Wall Clock Time and Write to File
#include <chrono>
#include <fstream>

// Base Components
#include "base.h"

// ME
#include "matrix.h"

// Shower
#include "shower.h"

// Analysis
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

void runGenerator(const int& N, const std::string& filename) {
  // ---------------------------------------------------------------------------
  // Give some information about the simulation

  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "|        GAPS: C++ Shower for Comparison        |" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Process: e+ e- --> q qbar" << std::endl;
  std::cout << "Number of Events: " << N << std::endl;
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // Inititalisation

  std::cout << "Initialising..." << std::endl;
  std::vector<Event> events(N);

  // ---------------------------------------------------------------------------
  // Matrix Element Generation

  std::cout << "Generating Matrix Elements (C++)..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  Matrix xs;

  for (int i = 0; i < N; i++) {
    xs.GenerateLOPoint(events[i]);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // Showering

  std::cout << "Showering Partons (C++)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  Shower sh;

  for (int i = 0; i < N; i++) {
    sh.Run(events[i]);
  }

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_sh = end - start;

  // ---------------------------------------------------------------------------
  // Analysis

  std::cout << "Analyzing Events (C++)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // Remove existing file
  std::remove(filename.c_str());

  Analysis an;

  // Analyze events (Including Validation of Colour and Momentum Conservation)
  for (int i = 0; i < N; i++) {
    an.Analyze(events[i]);
  }

  // Storage
  an.Finalize(filename);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_an = end - start;

  // ---------------------------------------------------------------------------
  // Results

  double diff = diff_me.count() + diff_sh.count() + diff_an.count();

  std::cout << "" << std::endl;
  std::cout << "EVENT GENERATION COMPLETE" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "ME Time: " << diff_me.count() << " s" << std::endl;
  std::cout << "Sh Time: " << diff_sh.count() << " s" << std::endl;
  std::cout << "An Time: " << diff_an.count() << " s" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Total Time: " << diff << " s" << std::endl;
  std::cout << "" << std::endl;

  // Open the file in append mode. This will create the file if it doesn't
  // exist.
  std::ofstream outfile("cpp-time.dat", std::ios_base::app);

  // Write diff_sh.count() to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // Close the file.
  outfile.close();

  std::cout << "Histograms written to " << filename << std::endl;
  std::cout << "Timing data written to cpp-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  int N = argc > 1 ? atoi(argv[1]) : 10000;
  runGenerator(N, "cpp.yoda");

  return 0;
}
// -----------------------------------------------------------------------------
