// to measure wall clock time and write to file
#include <chrono>
#include <fstream>

// base components
#include "base.cuh"

// matrix element
#include "matrix.cuh"

// parton shower
#include "shower.cuh"

// jet and event shape analysis
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

__global__ void set_seed_kernel(event* events, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx >= n) {
    return;
  }

  event& ev = events[idx];
  ev.set_seed(idx);
}

void run_generator(const int& n, const double& e, const std::string& filename) {
  // ---------------------------------------------------------------------------
  // Give some information about the simulation

  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "|      GAPS: a GPU-Amplified Parton Shower      |" << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Process: e+ e- --> q qbar" << std::endl;
  std::cout << "Number of Events: " << n << std::endl;
  std::cout << "Center of Mass Energy: " << e << " GeV" << std::endl;
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // initialisation

  std::cout << "Initialising..." << std::endl;
  thrust::device_vector<event> d_events(n);

  // set the seed
  set_seed_kernel<<<(n + 255) / 256, 256>>>(
      thrust::raw_pointer_cast(d_events.data()), n);

  // ---------------------------------------------------------------------------
  // me generation

  std::cout << "Generating matrix elements..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  calc_lome(d_events, e);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // do the showering

  std::cout << "Showering partons..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  run_shower(d_events);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_sh = end - start;

  // ---------------------------------------------------------------------------
  // analyze events

  std::cout << "Analysing events..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // analysis
  do_analysis(d_events, filename);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_an = end - start;

  // ---------------------------------------------------------------------------
  // empty the device vector (not neccessary, but good practice)

  d_events.clear();
  d_events.shrink_to_fit();

  /**
   * maybe in the future, to allow > 10^6 events, we can split the large number
   * into smaller batches. right now, we write events to file directly from the
   * do... functions, so the code is not ready for this.
   */

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
  std::cout << "Total Time: " << diff << " s" << std::endl;
  std::cout << "" << std::endl;

  // open the file in append mode. this will create the file if it doesn't
  // exist.
  std::ofstream outfile("gaps-time.dat", std::ios_base::app);

  // write diff_sh.count() to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // close the file.
  outfile.close();

  std::cout << "Histograms written to " << filename << std::endl;
  std::cout << "Timing data written to gaps-time.dat" << std::endl;
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
      std::string filename = "gaps-" + std::to_string(i) + ".yoda";
      run_generator(max_events, e, filename);
    }

    // run remainder
    if (n_remainder > 0) {
      std::string filename = "gaps-" + std::to_string(n_batches) + ".yoda";
      run_generator(n_remainder, e, filename);
    }
  } else {
    run_generator(n, e, "gaps.yoda");
  }

  return 0;
}
// -----------------------------------------------------------------------------
