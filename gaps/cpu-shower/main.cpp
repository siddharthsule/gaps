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

// interface for input
#include "interface.h"

/**
 * The Main Function
 * -----------------
 *
 * This file contains the main function for the CPU Shower program.
 * It is responsible for setting up the event generation, calling the
 * matrix element calculation, performing the parton showering, and
 * analyzing the final state particles.
 */

// -----------------------------------------------------------------------------

void run_generator(const params& p) {
  /**
   * @brief Run the event generator
   *
   * @param p The parameters for the run, passed from the main function
   */

  // ---------------------------------------------------------------------------
  // inititalisation

  std::cout << "Initialising..." << std::endl;
  std::vector<event> events(p.n_events);  // Create n events

  for (int i = 0; i < p.n_events; i++) {
    event& ev = events[i];
    ev.set_id(i + p.id_offset);  // Set the event id
    ev.set_seed(static_cast<unsigned long>(i + p.id_offset));  // Set the seed
    double dummy = ev.gen_random();  // Generate the first random number
  }

  std::cout << " - Using LHAPDF with ME2 PDF: " << p.me2pdf << std::endl;
  std::cout << " - Using LHAPDF with Shower PDF: " << p.showerpdf << std::endl;
  LHAPDF::setVerbosity(0);

  // Extra line to add space
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating Matrix Elements (CPU)..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  matrix me(p.process, p.nlo, p.root_s, p.asmz, p.me2pdf);

  for (int i = 0; i < p.n_events; i++) {
    me.run(events[i]);
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // showering

  std::chrono::duration<double> diff_sh(0.0);

  if (!p.no_shower) {
    std::cout << "Showering Partons (CPU)..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    shower sh(p.root_s, p.t_c, p.asmz, p.fixed_as, p.use_cmw, p.n_emissions_max,
              p.showerpdf);

    for (int i = 0; i < p.n_events; i++) {
      sh.run(events[i], p.nlo);
      std::cerr << "\rEvent " << i + 1 << " of " << p.n_events << std::flush;
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
  std::remove(p.storage_file.c_str());

  analysis an(p.process);

  // analyze events (including validation of colour and momentum conservation)
  for (int i = 0; i < p.n_events; i++) {
    an.validate(events[i]);
    an.analyze(events[i]);
  }
  std::cout << "" << std::endl;

  // storage
  an.finalize(p.storage_file);

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

  std::cout << "Histograms written to " << p.storage_file << std::endl;
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

  params run_params(argv);

  // run the generator
  run_generator(run_params);
  return 0;
}
// -----------------------------------------------------------------------------
