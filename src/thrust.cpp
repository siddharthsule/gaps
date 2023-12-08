#include "thrust.h"
#include <tuple>

Thrust::Thrust() : EShape() {}

std::tuple<double, std::vector<double>> Thrust::CalculateThrust(const Event& ev) {
  std::vector<std::vector<double>> moms = GetMomenta(ev);

  // Sum of all momenta
  double momsum = 0.0;
  for (auto& mom : moms) {
    momsum += Magnitude(mom);
  }

  double thrust = 0.0;
  std::vector<double> t_axis = {0.0, 0.0, 0.0};

  for (size_t k = 1; k < moms.size(); ++k) {
    for (size_t j = 0; j < k; ++j) {
      std::vector<double> tmp_axis = CrossProduct(moms[j], moms[k]);
      std::vector<double> p_thrust = {0.0, 0.0, 0.0};
      std::vector<std::vector<double>> p_combin;

      for (size_t i = 0; i < moms.size(); ++i) {
        if (i != j && i != k) {
          if (DotProduct(moms[i], tmp_axis) >= 0) {
            p_thrust = Add(p_thrust, moms[i]);
          } else {
            p_thrust = Subtract(p_thrust, moms[i]);
          }
        }
      }

      p_combin.push_back(Add(Add(p_thrust, moms[j]), moms[k]));
      p_combin.push_back(Subtract(Add(p_thrust, moms[j]), moms[k]));
      p_combin.push_back(Add(Subtract(p_thrust, moms[j]), moms[k]));
      p_combin.push_back(Subtract(Subtract(p_thrust, moms[j]), moms[k]));

      for (auto p : p_combin) {
        double temp = Magnitude(p);
        if (temp > thrust) {
          thrust = temp;
          t_axis = p;  // will unit-ify later
        }
      }
    }
  }

  thrust /= momsum;
  thrust = 1.0 - thrust;

  t_axis = Multiply(t_axis, 1.0 / Magnitude(t_axis));
  if (t_axis[2] < 0) {
    t_axis = Multiply(t_axis, -1.0);
  }

  return std::make_tuple(thrust, t_axis);
}

TAnalysis::TAnalysis() : thr(), thr_hist(100, 0., .5, "/GPU_EvGEN/tvalue\n"), wtot(0.0) {}

void TAnalysis::Analyze(const Event& ev) {
  
  std::tuple<double, std::vector<double>> thrust = thr.CalculateThrust(ev);

  double& tvalue = std::get<0>(thrust);
  if (tvalue < 1e-12) {
    tvalue = -5.0;
  }

  thr_hist.Fill(tvalue, ev.GetDxs());
  wtot += ev.GetDxs();
}

void TAnalysis::Finalize(const std::string& filename) {
  thr_hist.ScaleW(1.0 / wtot);
  thr_hist.Write(filename);

  Histo1D tzmHist = thr_hist;
  tzmHist.name = "/GPU_EvGEN/tzoomd";
  tzmHist.Write(filename);
}