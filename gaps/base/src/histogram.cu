#include "histogram.cuh"

// Libraries needed for file writing
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

// Text for the Yoda File
std::string to_string(histo1d h, std::string name) {
  std::stringstream ss;
  ss << "BEGIN YODA_HISTO1D " << name << "\n\n";
  ss << "Path=" << name << "\n\n";
  ss << "ScaledBy=" << h.scale << "\n";
  ss << "Title=\nType=Histo1D\n";
  ss << "# ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tnumEntries\n";
  ss << std::scientific << std::setprecision(6);
  ss << "Total"
     << "\t" << h.total.w << "\t" << h.total.w2 << "\t" << h.total.wx << "\t"
     << h.total.wx2 << "\t" << static_cast<int>(h.total.n) << "\n";
  ss << "Underflow"
     << "\t" << h.uflow.w << "\t" << h.uflow.w2 << "\t" << h.uflow.wx << "\t"
     << h.uflow.wx2 << "\t" << static_cast<int>(h.uflow.n) << "\n";
  ss << "Overflow"
     << "\t" << h.oflow.w << "\t" << h.oflow.w2 << "\t" << h.oflow.wx << "\t"
     << h.oflow.wx2 << "\t" << static_cast<int>(h.oflow.n) << "\n";
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

// write the yoda file
void write(histo1d h, std::string name, const std::string& filename) {
  std::ofstream file;
  file.open(filename, std::ios::out | std::ios::app);
  file << to_string(h, name);
  file.close();
}

std::string to_string(histo2d h, std::string name) {
  std::stringstream ss;
  ss << "BEGIN YODA_HISTO2D " << name << "\n\n";
  ss << "Path=" << name << "\n\n";
  ss << "ScaledBy=" << h.scale << "\n";
  ss << "Title=\nType=Histo2D\n";
  ss << "# "
        "ID\tID\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2\tsumwxy\tnumEntries"
        "\n";
  ss << std::scientific << std::setprecision(6);
  ss << "Total"
     << "\t" << h.total.w << "\t" << h.total.w2 << "\t" << h.total.wx << "\t"
     << h.total.wx2 << "\t" << h.total.wy << "\t" << h.total.wy2 << "\t"
     << h.total.wxy << "\t" << static_cast<int>(h.total.n) << "\n";
  ss << "# "
        "xlow\txhigh\tylow\tyhigh\tsumw\tsumw2\tsumwx\tsumwx2\tsumwy\tsumwy2\ts"
        "umwxy\tnumEntries\n";
  for (size_t i = 0; i < n_bins2d; ++i) {
    for (size_t j = 0; j < n_bins2d; ++j) {
      ss << std::scientific << std::setprecision(6);
      ss << h.bins[i][j].xmin << "\t" << h.bins[i][j].xmax << "\t"
         << h.bins[i][j].ymin << "\t" << h.bins[i][j].ymax << "\t"
         << h.bins[i][j].w << "\t" << h.bins[i][j].w2 << "\t" << h.bins[i][j].wx
         << "\t" << h.bins[i][j].wx2 << "\t" << h.bins[i][j].wy << "\t"
         << h.bins[i][j].wy2 << "\t" << h.bins[i][j].wxy << "\t"
         << static_cast<int>(h.bins[i][j].n) << "\n";
    }
  }
  ss << "END YODA_HISTO2D\n\n";
  return ss.str();
}

void write(histo2d h, std::string name, const std::string& filename) {
  std::ofstream file;
  file.open(filename, std::ios::out | std::ios::app);
  file << to_string(h, name);
  file.close();
}