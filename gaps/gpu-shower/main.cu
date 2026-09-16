// to measure wall clock time and write to file
#include <chrono>
#include <fstream>

// base components
#include "base.cuh"

// matrix element
#include "matrix.cuh"

// parton shower
#include "shower.cuh"

// hadronisation
#include "hadronisation.cuh"

// decays
#include "hadronic_decays.cuh"

// jet and event shape analysis
#include "observables.cuh"

// interface for input
#include "interface.cuh"

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
// kernel to set the seed for the random number generator

__global__ void set_seed_kernel(event* events, int id_offset, int n) {
  /**
   * @brief Set the seed for the random number generator
   *
   * @param events array of event records
   * @param id_offset the offset for the event id
   * @param n number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------

  event& ev = events[idx];
  ev.set_id(idx + id_offset);
  ev.set_seed(static_cast<unsigned long>(idx + id_offset));
  double dummy = ev.gen_random();
}

void run_generator(const params& p) {
  /**
   * @brief Run the event generator
   *
   * @param p The parameters for the run, passed from the main function
   */
  // ---------------------------------------------------------------------------
  // initialisation

  std::cout << "Initialising..." << std::endl;

  // create the events
  thrust::device_vector<event> dv_events(p.n_events);
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // Threads-per-Block and Blocks-per-Grid
  int blocks =
      static_cast<int>(std::ceil(static_cast<double>(n_events) / p.threads));
  std::cout << " - Using " << blocks << " blocks and " << p.threads
            << " threads per block" << std::endl;

  // set the seed
  set_seed_kernel<<<blocks, p.threads>>>(d_events, p.id_offset, p.n_events);

  // Output LHAPDF settings
  std::cout << " - Using LHAPDF with ME2 PDF: " << p.me2pdf << std::endl;
  std::cout << " - Using LHAPDF with Shower PDF: " << p.showerpdf << std::endl;
  LHAPDF::setVerbosity(0);

  // Extra line to add space
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating matrix elements..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  // Calculate the leading order cross section and kinematics
  run_matrix(dv_events, p, blocks);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // do the showering

  std::chrono::duration<double> diff_sh(0.0);

  if (!p.no_shower) {
    std::cout << "Showering partons..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    run_shower(dv_events, p, blocks);

    end = std::chrono::high_resolution_clock::now();
    diff_sh = end - start;
  } else {
    std::cout << "Skipping shower section (noshower enabled)..." << std::endl;
  }

  // ---------------------------------------------------------------------------
  // hadronisation

  std::chrono::duration<double> diff_had(0.0);

  if (p.hadronise) {
    // Cluster lists, one per event. Allocated before the timer so
    // their initialisation is not counted as hadronisation time.
    thrust::device_vector<cluster_list> dv_cls(p.n_events);

    std::cout << "Hadronising clusters..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    run_hadronisation(dv_events, dv_cls, p, blocks);

    end = std::chrono::high_resolution_clock::now();
    diff_had = end - start;
  } else {
    std::cout << "Skipping hadronisation section (hadronise disabled)..."
              << std::endl;
  }

  // ---------------------------------------------------------------------------
  // hadronic decays

  std::chrono::duration<double> diff_dec(0.0);

  if (p.hadronise) {
    std::cout << "Decaying Hadrons..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    run_decays(dv_events, p, blocks);

    end = std::chrono::high_resolution_clock::now();
    diff_dec = end - start;
  } else {
    std::cout << "Skipping decay section (hadronise disabled)..." << std::endl;
  }

  // ---------------------------------------------------------------------------
  // analyze events

  std::chrono::duration<double> diff_an(0.0);

  if (!p.skip_analysis) {
    std::cout << "Analysing events..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    do_analysis(dv_events, p, blocks);

    end = std::chrono::high_resolution_clock::now();
    diff_an = end - start;
  } else {
    std::cout << "Skipping analysis section (skip_analysis enabled)..."
              << std::endl;
  }

  // ---------------------------------------------------------------------------
  // Additionally Try Moving Events to Host

  // std::cout << "Additional Test: Moving Events to Host..." << std::endl;
  // start = std::chrono::high_resolution_clock::now();

  // thrust::host_vector<event> h_events = dv_events;

  // for (int i = 0; i < n; i++) {
  //   if (!h_events[i].get_validity()) {
  //     std::cout << "Invalid Event: " << i << std::endl;
  //     h_events[i].print_info();
  //   }
  // }

  // h_events[0].print_info();  // print the first event

  // end = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> diff_host = end - start;

  // ---------------------------------------------------------------------------
  // results

  double diff = diff_me.count() + diff_sh.count() + diff_had.count() +
                diff_dec.count() + diff_an.count();

  std::cout << "" << std::endl;
  std::cout << "EVENT GENERATION COMPLETE" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "ME Time: " << diff_me.count() << " s" << std::endl;
  std::cout << "Sh Time: " << diff_sh.count() << " s" << std::endl;
  std::cout << "Hd Time: " << diff_had.count() << " s" << std::endl;
  std::cout << "Dc Time: " << diff_dec.count() << " s" << std::endl;
  std::cout << "An Time: " << diff_an.count() << " s" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Total Time: " << diff << " s" << std::endl;
  std::cout << "" << std::endl;
  // std::cout << "Moving Events to Host: " << diff_host.count() << " s"
  //           << std::endl;

  // open the file in append mode. this will create the file if it doesn't
  // exist.
  std::ofstream outfile("gpu-time.dat", std::ios_base::app);

  // write the timings to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_had.count() << ", " << diff_dec.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // close the file.
  outfile.close();

  if (!p.skip_analysis) {
    std::cout << "Histograms written to " << p.storage_file << std::endl;
  }
  std::cout << "Timing data written to gpu-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  /**
   * @brief Main function to run the GPU Shower
   *
   * All validation is done in the Python interface.
   */

  params run_params(argv);

  // more than max_events must be split into batches by the Python interface
  if (run_params.n_events > max_events) {
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    return 1;
  }

  // run the generator
  run_generator(run_params);
  return 0;
}
// -----------------------------------------------------------------------------
