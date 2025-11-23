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
 * The Main Function
 * -----------------
 *
 * This file contains the main function for the GPU Shower program.
 * It is responsible for setting up the event generation, calling the
 * matrix element calculation, performing the parton showering, and
 * analyzing the final state particles.
 */

// -----------------------------------------------------------------------------

void run_generator(int process, bool nlo, double root_s, double asmz,
                   bool fixed_as, bool no_shower, double t_c,
                   int n_emissions_max, int n, int id_offset,
                   std::string filename) {
  /**
   * @brief Run the event generator
   *
   * @param process: Process ID
   * @param nlo: NLO or LO
   * @param root_s: Center of mass energy
   * @param asmz: Strong coupling constant at Z mass
   * @param fixed_as: Whether to use fixed strong coupling
   * @param no_shower: Whether to skip the shower section
   * @param t_c: Shower cutoff in GeV
   * @param n_emissions_max: Maximum number of emissions
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

  std::cout << " - Using LHAPDF with CT14lo set" << std::endl;
  LHAPDF::setVerbosity(0);

  // Extra line to add space
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating Matrix Elements (CPU)..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  matrix me(process, nlo, root_s, asmz);

  for (int i = 0; i < n; i++) {
    me.run(events[i]);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // showering

  std::chrono::duration<double> diff_sh(0.0);

  if (!no_shower) {
    std::cout << "Showering Partons (CPU)..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    shower sh(root_s, t_c, asmz, fixed_as, n_emissions_max);

    for (int i = 0; i < n; i++) {
      sh.run(events[i], nlo);
      std::cerr << "\rEvent " << i + 1 << " of " << n << std::flush;
    }
    std::cout << "" << std::endl;

    end = std::chrono::high_resolution_clock::now();
    diff_sh = end - start;
  } else {
    std::cout << "Skipping shower section (noshower enabled)..." << std::endl;
  }

  // ---------------------------------------------------------------------------
  // analysis

  std::cout << "Analysing Events (CPU)..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // remove existing file
  std::remove(filename.c_str());

  analysis an(process);

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

  // 1 - 1 for LEP or 2 for LHC
  int process = atoi(argv[1]);

  // 2 - 0 if no NLO, 1 if NLO
  bool nlo = atoi(argv[2]);

  // 3 - The root s value
  double root_s = atof(argv[3]);

  // 4 - The strong coupling constant at the Z mass
  double asmz = atof(argv[4]);

  // 5 - If fixas is set, use the fixed asmz value
  bool fixed_as = atoi(argv[5]);

  // 6 - If noshower is set, skip the shower section
  bool no_shower = atoi(argv[6]);

  // 7 - The Shower Cutoff in GeV
  double t_c = atof(argv[7]);

  // 8 - The maximum number of emissions
  int n_emissions_max = atoi(argv[8]);

  // 9 - The Number of events
  int n_events = atoi(argv[9]);

  // 10 - The Event Number Offset
  int id_offset = atoi(argv[10]);

  // 11 - Storage file name
  std::string storage_file = argv[11];

  // if more than max_events, run in batches
  if (n_events > max_events) {
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    return 1;
  }

  // run the generator
  run_generator(process, nlo, root_s, asmz, fixed_as, no_shower, t_c,
                n_emissions_max, n_events, id_offset, storage_file);
  return 0;
}
// -----------------------------------------------------------------------------
