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

void run_generator(bool nlo, double root_s, double asmz, double t_c,
                   int n_emissions_max, int n, int id_offset,
                   std::string filename, const int& threads) {
  /**
   * @brief Run the event generator
   *
   * @param nlo: NLO or LO
   * @param root_s: Center of mass energy
   * @param n: Number of events to generate
   * @param id_offset: Offset for event IDs
   * @param filename: Name of the file to store the histograms
   * @param threads: Number of threads per block
   */
  // ---------------------------------------------------------------------------
  // initialisation

  std::cout << "Initialising..." << std::endl;

  // create the events
  thrust::device_vector<event> dv_events(n);
  event* d_events = thrust::raw_pointer_cast(dv_events.data());
  int n_events = dv_events.size();

  // Threads-per-Block and Blocks-per-Grid
  // int threads = 256;
  int blocks =
      static_cast<int>(std::ceil(static_cast<double>(n_events) / threads));
  std::cout << " - Using " << blocks << " blocks and " << threads
            << " threads per block." << std::endl;

  // set the seed
  set_seed_kernel<<<blocks, threads>>>(d_events, id_offset, n);

  // Extra line to add space
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating matrix elements..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  // Calculate the leading order cross section and kinematics
  calc_lome(dv_events, nlo, root_s, asmz, blocks, threads);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  // ---------------------------------------------------------------------------
  // do the showering

  std::cout << "Showering partons..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  bool do_partition = true;
  run_shower(dv_events, root_s, nlo, do_partition, t_c, asmz, n_emissions_max,
             blocks, threads);

  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_sh = end - start;

  // ---------------------------------------------------------------------------
  // analyze events

  std::cout << "Analysing events..." << std::endl;
  start = std::chrono::high_resolution_clock::now();

  // analysis
  do_analysis(dv_events, filename, blocks, threads);

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

  std::cout << "Histograms written to " << filename << std::endl;
  std::cout << "Timing data written to gpu-time.dat" << std::endl;
  std::cout << "------------------------------------------------" << std::endl;

  // ---------------------------------------------------------------------------
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

  // Threads per block
  int threads = atoi(argv[9]);

  // if more than max_events, run in batches
  if (n_events > max_events) {
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    return 1;
  }

  // run the generator
  run_generator(nlo, root_s, asmz, t_c, n_emissions_max, n_events, id_offset,
                storage_file, threads);
  return 0;
}
// -----------------------------------------------------------------------------
