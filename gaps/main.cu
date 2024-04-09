// To Measure Wall Clock Time and Write to File
#include <chrono>
#include <fstream>

// Base Components
#include "base.cuh"

// Matrix Element
#include "matrix.cuh"

// Parton Shower
#include "shower.cuh"

// Jet and Event Shape Analysis
#include "observables.cuh"

/**
 * GAPS: a GPU-Amplified Parton Shower
 * -----------------------------------
 *
 * This program is a simple event generator for e+e- -> partons. It is designed
 * to be a simple example of how to use the GPU to calculate matrix elements and
 * perform parton showering. The program is designed to be a proof of concept as
 * well as a intuitive example of how to use the GPU for event generation.
 *
 * This program is based on S. Höche's "Introduction to Parton Showers" Python
 * tutorial[1]. This program calculates ee -> qq and then showers the partons.
 *
 * [1] https://arxiv.org/abs/1411.4085 and MCNET-CTEQ 2021 Tutorial
 */

// -----------------------------------------------------------------------------

void runGenerator(const int& N, const double& E, const std::string& filename) {
  // ---------------------------------------------------------------------------
  // Give some information about the simulation

  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "|      GAPS: a GPU-Amplified Parton Shower      |" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Process: e+ e- --> q qbar" << std::endl;
  std::cout << "Number of Events: " << N << std::endl;
  std::cout << "Center of Mass Energy: " << E << " GeV" << std::endl;
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // Initialisation

  std::cout << "Initialising..." << std::endl;
  thrust::device_vector<Event> d_events(N);

  // ---------------------------------------------------------------------------
  // ME Generation

  std::cout << "Generating Matrix Elements..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  calcLOME(d_events, E);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // Do the Showering

  std::cout << "Showering Partons..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  runShower(d_events);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_sh = end - start;

  // ---------------------------------------------------------------------------
  // Analyze Events

  std::cout << "Analyzing Events..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // Analysis
  doAnalysis(d_events, filename);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_an = end - start;

  // ---------------------------------------------------------------------------
  // Empty the device vector (Not neccessary, but good practice)

  d_events.clear();
  d_events.shrink_to_fit();

  /**
   * Maybe in the future, to allow > 10^6 events, we can split the large number
   * into smaller batches. Right now, we write events to file directly from the
   * do... functions, so the code is not ready for this.
   */

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
  std::ofstream outfile("gaps-time.dat", std::ios_base::app);

  // Write diff_sh.count() to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // Close the file.
  outfile.close();

  std::cout << "Histograms written to " << filename << std::endl;
  std::cout << "Timing data written to gaps-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  // Import Settings
  int N = argc > 1 ? atoi(argv[1]) : 10000;
  double E = argc > 2 ? atof(argv[2]) : 91.2;

  // If more than maxEvents, run in batches
  if (N > maxEvents) {
    std::cout << "-------------------------------------------------"
              << std::endl;
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    std::cout << "Running in batches..." << std::endl;
    std::cout << "Please use rivet-merge to combine runs" << std::endl;

    // Split into batches
    int nBatches = N / maxEvents;
    int nRemainder = N % maxEvents;
    std::cout << "Number of Batches: " << nBatches << std::endl;
    std::cout << "Size of Remainder: " << nRemainder << std::endl;

    // Run in batches
    for (int i = 0; i < nBatches; i++) {
      std::string filename = "gaps-" + std::to_string(i) + ".yoda";
      runGenerator(maxEvents, E, filename);
    }

    // Run remainder
    if (nRemainder > 0) {
      std::string filename = "gaps-" + std::to_string(nBatches) + ".yoda";
      runGenerator(nRemainder, E, filename);
    }
  } else {
    runGenerator(N, E, "gaps.yoda");
  }

  return 0;
}
// -----------------------------------------------------------------------------
