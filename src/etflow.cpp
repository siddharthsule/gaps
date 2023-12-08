#include "etflow.h"
#include <tuple>

ETFlow::ETFlow() {}

std::vector<double> ETFlow::CalculateETFlow(const Event& ev) {
  
  std::vector<std::vector<double>> moms = GetMomenta(ev);

  std::tuple<double, std::vector<double>> thrust = thr.CalculateThrust(ev);
  std::vector<double> t_axis = std::get<1>(thrust);

  std::vector<double> etflow = {0.0, 0.0, 0.0};

  for (const auto& mo : moms) {
    double mo_para = DotProduct(mo, t_axis);
    double mo_perp = Magnitude(Subtract(mo, Multiply(t_axis, mo_para)));
    double enrg = Magnitude(mo);

    try {
      double y = 0.5 * std::log((enrg - mo_para) / (enrg + mo_para));
      y = std::fabs(y);

      if (y < 0.5) {
        etflow[0] += mo_perp;
      }
      if (y < 1.0) {
        etflow[1] += mo_perp;
      }
      if (y < 1.5) {
        etflow[2] += mo_perp;
      }
    } catch (const std::exception& e) {
      // y -> +- inf, ie negligible Et
      // std::cout << "parallel to thrust" << std::endl;
      continue;
    }
  }

  return etflow;
}

ETAnalysis::ETAnalysis() : etflow(), wtot(0.0) {
  etflow_hists.push_back(Histo1D(100, 0.0, 80.0, "/GPU_EvGEN/et0p5\n"));
  etflow_hists.push_back(Histo1D(100, 0.0, 80.0, "/GPU_EvGEN/et1p0\n"));
  etflow_hists.push_back(Histo1D(100, 0.0, 80.0, "/GPU_EvGEN/et1p5\n"));
}

void ETAnalysis::Analyze(const Event& ev) {
  std::vector<double> etdata = etflow.CalculateETFlow(ev);

  etflow_hists[0].Fill(etdata[0], ev.GetDxs());
  etflow_hists[1].Fill(etdata[1], ev.GetDxs());
  etflow_hists[2].Fill(etdata[2], ev.GetDxs());

  wtot += ev.GetDxs();
}

void ETAnalysis::Finalize(const std::string& filename) {
  for (auto& hist : etflow_hists) {
    hist.ScaleW(1.0 / wtot);
    hist.Write(filename);
  }
}


