#ifndef interface_h_
#define interface_h_

#include <string>

class params {
  /**
   * @class params
   * @brief The parameters for the codebase, kept in a class for easy
   * transfer between different functions and classes.
   */

 public:
  int process;               // 1 = ee2qq, 11 = pp2Z, 12 = pp2yy
  bool nlo;                  // Whether to run in NLO mode
  double root_s;             // Center of mass energy
  double asmz;               // Strong coupling constant at Z mass
  bool fixed_as;             // Whether to use fixed strong coupling
  bool use_cmw;              // Whether to apply CMW scheme conversion
  bool no_shower;            // Whether to skip the shower section
  double t_c;                // Shower cutoff in GeV
  int n_emissions_max;       // Maximum number of emissions
  std::string me2pdf;        // PDF set name for matrix element calculation
  std::string showerpdf;     // PDF set name for parton shower
  int n_events;              // Number of events to generate
  int id_offset;             // Offset for event IDs
  std::string storage_file;  // Name of the file to store the histograms

  // Constructor
  params(char* argv[]) {
    process = atoi(argv[1]);
    nlo = atoi(argv[2]);
    root_s = atof(argv[3]);
    asmz = atof(argv[4]);
    fixed_as = atoi(argv[5]);
    use_cmw = atoi(argv[6]);
    no_shower = atoi(argv[7]);
    t_c = atof(argv[8]);
    n_emissions_max = atoi(argv[9]);
    me2pdf = argv[10];
    showerpdf = argv[11];
    n_events = atoi(argv[12]);
    id_offset = atoi(argv[13]);
    storage_file = argv[14];
  }
};

#endif  // interface_h_