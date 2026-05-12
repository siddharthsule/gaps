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
  calc_lome(dv_events, p.process, p.nlo, p.root_s, p.asmz, blocks, p.threads,
            p.me2pdf);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  thrust::host_vector<event> h_events = dv_events;

  // ---------------------------------------------------------------------------
  // do the showering

  std::chrono::duration<double> diff_sh(0.0);

  if (!p.no_shower) {
    std::cout << "Showering partons..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    run_shower(dv_events, p.root_s, p.nlo, p.do_partitioning, p.t_c, p.asmz,
               p.fixed_as, p.use_cmw, p.n_emissions_max, blocks, p.threads,
               p.showerpdf);

    end = std::chrono::high_resolution_clock::now();
    diff_sh = end - start;
  } else {
    std::cout << "Skipping shower section (noshower enabled)..." << std::endl;
  }

  // ---------------------------------------------------------------------------
  // analyze events

  std::cout << "Analysing events..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // analysis
  do_analysis(dv_events, p.storage_file, p.process, blocks, p.threads);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_an = end - start;

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
  // std::cout << "Moving Events to Host: " << diff_host.count() << " s"
  //           << std::endl;

  // open the file in append mode. this will create the file if it doesn't
  // exist.
  std::ofstream outfile("gpu-time.dat", std::ios_base::app);

  // write diff_sh.count() to the file.
  outfile << diff_me.count() << ", " << diff_sh.count() << ", "
          << diff_an.count() << ", " << diff << std::endl;

  // close the file.
  outfile.close();

  std::cout << "Histograms written to " << p.storage_file << std::endl;
  std::cout << "Timing data written to gpu-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
}
// -----------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  /**
   * @brief Main function to run the GPU Shower
   *
   * All Validation is done in the Python Interface, so here is just the
   * main function to run the generator. We simply add one check for the
   * number of events.
   */

  params run_params(argv);

  // if more than max_events, run in batches
  if (run_params.n_events > max_events) {
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    return 1;
  }

  // run the generator
  run_generator(run_params);
  return 0;
}
// -----------------------------------------------------------------------------
