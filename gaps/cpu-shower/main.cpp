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
 * GAPS: CPU Shower for Comparison
 * ------------------------------------
 *
 * This program is a translation of S. Höche's "Introduction to Parton Showers"
 * Python tutorial[1], with added functionality for parallelisation, a Event
 * class and event shape analyses.
 *
 * The purpose of this program is to compare the performance of the CPU and
 * GPU versions of the shower, and to compare the performance of the CPU with
 * parallelisation and GPU.
 *
 * [1] https://arxiv.org/abs/1411.4085 and MCNET-CTEQ 2021 Tutorial
 */

// -----------------------------------------------------------------------------

void run_generator(bool nlo, double root_s, double asmz, double t_c,
                   int n_emissions_max, int n, int id_offset,
                   std::string filename) {
  /**
   * @brief Run the event generator
   *
   * @param nlo: NLO or LO
   * @param root_s: Center of mass energy
   * @param n: Number of events to generate
   * @param id_offset: Offset for event IDs
   * @param filename: Name of the file to store the histograms
   */

  // ---------------------------------------------------------------------------
  // inititalisation

  std::cout << "Initialising..." << std::endl;
  std::vector<event> events(n);  // Create n events

  for (int i = 0; i < n; i++) {
    event& ev = events[i];
    ev.set_id(i + id_offset);                                // Set the event id
    ev.set_seed(static_cast<unsigned long>(i + id_offset));  // Set the seed
    double dummy = ev.gen_random();  // Generate the first random number
  }

  // Extra line to add space
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating Matrix Elements (CPU)..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  matrix me(nlo, root_s, asmz);

  for (int i = 0; i < n; i++) {
    me.run(events[i]);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // showering

  std::cout << "Showering Partons (CPU)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  shower sh(root_s, t_c, asmz, n_emissions_max);

  for (int i = 0; i < n; i++) {
    sh.run(events[i], nlo);
    std::cerr << "\rEvent " << i + 1 << " of " << n << std::flush;
  }
  std::cout << "" << std::endl;

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_sh = end - start;

  // ---------------------------------------------------------------------------
  // analysis

  std::cout << "Analysing Events (CPU)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // remove existing file
  std::remove(filename.c_str());

  analysis an;

  // analyze events (including validation of colour and momentum conservation)
  for (int i = 0; i < n; i++) {
    an.validate(events[i]);
    an.analyze(events[i]);
  }
  std::cout << "" << std::endl;

  // storage
  an.finalize(filename);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_an = end - start;

  // ---------------------------------------------------------------------------
  // results

  // events[0].print_info();  // print the first event

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
  std::ofstream outfile("cpu-time.dat", std::ios_base::app);

  // write diff_sh.count() to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // close the file.
  outfile.close();

  std::cout << "Histograms written to " << filename << std::endl;
  std::cout << "Timing data written to cpu-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  /**
   * @brief Main function to run the CPU Shower
   *
   * All Validation is done in the Python Interface, so here is just the
   * main function to run the generator. We simply add one check for the
   * number of events.
   */

  // NLO or LO? O = false, 1 = true
  bool nlo = atoi(argv[1]);

  // Energy or Root_S
  double root_s = atof(argv[2]);

  // Alpha(s) at Z mass
  double asmz = atof(argv[3]);

  // Shower Cutoff in GeV
  double t_c = atof(argv[4]);

  // Max number of emission
  int n_emissions_max = atoi(argv[5]);

  // Number of events
  int n_events = atoi(argv[6]);

  // Event ID Offset
  int id_offset = atoi(argv[7]);

  // Storage file name
  std::string storage_file = argv[8];

  // if more than max_events, run in batches
  if (n_events > max_events) {
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    return 1;
  }

  // run the generator
  run_generator(nlo, root_s, asmz, t_c, n_emissions_max, n_events, id_offset,
                storage_file);
  return 0;
}
// -----------------------------------------------------------------------------
