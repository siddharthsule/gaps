#include "histogram.cuh"

// Libraries needed for file writing
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

std::string to_string(histo1d h) {
  /**
   * @brief convert the histogram data to a string
   *
   * @param h the histogram to convert
   * @return std::string: the histogram data as a string
   */

  // Convert the char array name to a std::string
  std::string name(h.name);

  std::stringstream ss;
  ss << "BEGIN YODA_HISTO1D " << name << "\n\n";
  ss << "Path=" << name << "\n\n";
  ss << "ScaledBy=" << h.scale << "\n";
  ss << "Title=\nType=Histo1D\n";
  ss << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
  ss << std::scientific << std::setprecision(6);
  ss << "Total\tTotal" << "\t" << h.total.w << "\t" << h.total.w2 << "\t"
     << h.total.wx << "\t" << h.total.wx2 << "\t" << static_cast<int>(h.total.n)
     << "\n";
  ss << "Underflow\tUnderflow" << "\t" << h.uflow.w << "\t" << h.uflow.w2
     << "\t" << h.uflow.wx << "\t" << h.uflow.wx2 << "\t"
     << static_cast<int>(h.uflow.n) << "\n";
  ss << "Overflow\tOverflow" << "\t" << h.oflow.w << "\t" << h.oflow.w2 << "\t"
     << h.oflow.wx << "\t" << h.oflow.wx2 << "\t" << static_cast<int>(h.oflow.n)
     << "\n";
  ss << "# xlow\txhigh\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
  for (size_t i = 0; i < n_bins; ++i) {
    ss << std::scientific << std::setprecision(6);
    ss << h.bins[i].xmin << "\t" << h.bins[i].xmax << "\t" << h.bins[i].w
       << "\t" << h.bins[i].w2 << "\t" << h.bins[i].wx << "\t" << h.bins[i].wx2
       << "\t" << static_cast<int>(h.bins[i].n) << "\n";
  }
  ss << "END YODA_HISTO1D\n\n";
  return ss.str();
}

void write(histo1d h, const std::string& filename) {
  /**
   * @brief write the histogram data to a yoda file
   *
   * @param h the histogram to write
   * @param filename the name of the file to write the histogram data to
   */

  std::ofstream file;
  file.open(filename, std::ios::out | std::ios::app);
  file << to_string(h);
  file.close();
}