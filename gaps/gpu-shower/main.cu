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

void run_generator(int process, bool nlo, double root_s, double asmz,
                   bool fixed_as, bool no_shower, double t_c,
                   int n_emissions_max, int n, int id_offset,
                   std::string filename, bool do_partitioning, int threads) {
  /**
   * @brief Run the event generator
   *
   * @param process the process to generate (1 for LEP, 2 for LHC)
   * @param nlo whether to use NLO matching
   * @param root_s the center-of-mass energy
   * @param asmz the strong coupling constant at the Z mass
   * @param fixed_as whether to use a fixed strong coupling constant
   * @param no_shower whether to skip the shower section
   * @param t_c the shower cutoff in GeV
   * @param n_emissions_max the maximum number of emissions
   * @param n the number of events to generate
   * @param id_offset the offset for the event id
   * @param filename the name of the file to store the histograms
   * @param do_partitioning whether to use event partitioning
   * @param threads the number of threads per block
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
            << " threads per block" << std::endl;

  // set the seed
  set_seed_kernel<<<blocks, threads>>>(d_events, id_offset, n);

  // Output LHAPDF settings
  std::cout << " - Using LHAPDF with CT14lo set" << std::endl;
  LHAPDF::setVerbosity(0);

  // Extra line to add space
  std::cout << "" << std::endl;

  // ---------------------------------------------------------------------------
  // matrix element generation

  std::cout << "Generating matrix elements..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  // Calculate the leading order cross section and kinematics
  calc_lome(dv_events, process, nlo, root_s, asmz, blocks, threads);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff_me = end - start;

  thrust::host_vector<event> h_events = dv_events;

  // ---------------------------------------------------------------------------
  // do the showering

  std::chrono::duration<double> diff_sh(0.0);

  if (!no_shower) {
    std::cout << "Showering partons..." << std::endl;
    start = std::chrono::high_resolution_clock::now();

    run_shower(dv_events, root_s, nlo, do_partitioning, t_c, asmz, fixed_as,
               n_emissions_max, blocks, threads);

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
  do_analysis(dv_events, filename, process, blocks, threads);

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

  // 12 - Do partitioning
  bool do_partitioning = atoi(argv[12]);

  // 13 - Threads per block
  int threads = atoi(argv[13]);

  // if more than max_events, run in batches
  if (n_events > max_events) {
    std::cout << "More Events than GPU Can Handle at Once!" << std::endl;
    return 1;
  }

  // run the generator
  run_generator(process, nlo, root_s, asmz, fixed_as, no_shower, t_c,
                n_emissions_max, n_events, id_offset, storage_file,
                do_partitioning, threads);
  return 0;
}
// -----------------------------------------------------------------------------
