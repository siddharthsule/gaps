#ifndef PARTON_H_
#define PARTON_H_

// Partons have Vec4 Momentum, Vec4 #includes Base
#include "vec4.h"

class Parton {
 public:
  // Constructor
  // Used at Multiple Points [HOST + DEVICE]
  Parton(int pid = 0, Vec4 momentum = Vec4(), int col = 0, int anticol = 0)
      : pid(pid), mom(momentum), col(col), anticol(anticol) {}

  // Getters and Setters
  int GetPid() const { return pid; }
  Vec4 GetMom() const { return mom; }
  int GetCol() const { return col; }
  int GetAntiCol() const { return anticol; }

  void SetPid(int pid) { this->pid = pid; }
  void SetMom(Vec4 mom) { this->mom = mom; }
  void SetCol(int col) { this->col = col; }
  void SetAntiCol(int anticol) { this->anticol = anticol; }

  // Boolean - If two partons are in a Colour Connected Dipole
  bool IsColorConnected(Parton p) {
    return (col > 0 && col == p.anticol) || (anticol > 0 && anticol == p.col);
  }

 private:
  int pid;
  Vec4 mom;
  int col;
  int anticol;
};

// Event Class
// Built to contain the partons and the dxs as one accesible object
// In future, can use to store thrust, log10y23 to parallelise those
class Event {
 private:
  // Temporary Solution - Allows 30 Partons
  // Mainly done to avoid using std::vector for GPU Showering
  // But of course, this is not a good solution
  Parton partons[maxPartons];

  // ME Params -----------------------------------------------------------------

  double dxs;  // Differential Cross Section
  int nHard;   // Number of Hard Partons

  // Shower Params -------------------------------------------------------------

  int nEmission;                     // Number of Emissions
  double showerT, showerZ, showerY;  // Evolution and Splitting Variables
  int showerC;                       // Colour Counter

  // Selecting Winner Emission - Defaults Values which represent no winner
  int winSF = 16;
  int winDipole[2] = {-1, -1};
  double winParams[2] = {0.0, 0.0};

  bool endShower = false;  // Shower End Flag - used if T < 1 GeV

  // Analysis Variables --------------------------------------------------------

  // Event Validity - Momentum and Colour Conservation
  bool validity = true;

  // Jet Rates using the Durham Algorithm
  double y23 = -50.0, y34 = -50.0, y45 = -50.0, y56 = -50.0;

  // Event Shape Variables - Thrust, Jet Masses and Broadenings
  double thr = -50.0, hjm = -50.0, ljm = -50.0, wjb = -50.0, njb = -50.0;
  Vec4 t_axis = Vec4();

 public:
  // Constructor ---------------------------------------------------------------

  // Empty, so that we can build our ME, PS onto it
  Event() {}

  // Getters -------------------------------------------------------------------

  // Access Partons in the Event
  Parton GetParton(int i) const { return partons[i]; }
  int GetSize() const { return nHard + nEmission; }
  int GetHard() const { return nHard; }
  int GetEmissions() const { return nEmission; }
  int GetPartonSize() const { return nHard + nEmission - 2; }

  // Get Differential Cross Section
  double GetDxs() const { return dxs; }


  // Get Shower Params
  double GetShowerT() const { return showerT; }
  double GetShowerY() const { return showerY; }
  double GetShowerZ() const { return showerZ; }
  int GetShowerC() const { return showerC; }

  // Get Winner Emission
  int GetWinSF() const { return winSF; }
  int GetWinDipole(int i) const { return winDipole[i]; }
  double GetWinParam(int i) const { return winParams[i]; }

  // Get Analysis Variables
  bool GetValidity() const { return validity; }

  double GetY23() const { return y23; }
  double GetY34() const { return y34; }
  double GetY45() const { return y45; }
  double GetY56() const { return y56; }
  double GetThr() const { return thr; }
  double GetHJM() const { return hjm; }
  double GetLJM() const { return ljm; }
  double GetWJB() const { return wjb; }
  double GetNJB() const { return njb; }

  Vec4 GetTAxis() const { return t_axis; }

  // Setters -------------------------------------------------------------------

  // Add / Replace Parton
  void SetParton(int i, Parton parton) { partons[i] = parton; }

  // Not used in ME [HOST]
  void SetPartonPid(int i, int pid) { partons[i].SetPid(pid); }
  void SetPartonMom(int i, Vec4 mom) { partons[i].SetMom(mom); }
  void SetPartonCol(int i, int col) { partons[i].SetCol(col); }
  void SetPartonAntiCol(int i, int anticol) { partons[i].SetAntiCol(anticol); }

  // Set Differential Cross Section and nHard
  void SetDxs(double dxs) { this->dxs = dxs; }
  void SetHard(int nHard) { this->nHard = nHard; }

  // Adjust and Increment Number of Emissions
  void SetEmissions(int nEmission) { this->nEmission = nEmission; }
  void IncrementEmissions() { nEmission++; }

  // Set Shower Params
  void SetShowerT(double showerT) { this->showerT = showerT; }
  void SetShowerY(double showerY) { this->showerY = showerY; }
  void SetShowerZ(double showerZ) { this->showerZ = showerZ; }

  void SetShowerC(int showerC) { this->showerC = showerC; }
  void IncrementShowerC() { showerC++; }

  // Set Winner Emission
  void SetWinSF(int winSF) { this->winSF = winSF; }
  void SetWinDipole(int i, int winDipole) { this->winDipole[i] = winDipole; }
  void SetWinParam(int i, double winParams) { this->winParams[i] = winParams; }

  // Set Analysis Variables
  void SetValidity(bool validity) { this->validity = validity; }

  void SetY23(double y23) { this->y23 = y23; }
  void SetY34(double y34) { this->y34 = y34; }
  void SetY45(double y45) { this->y45 = y45; }
  void SetY56(double y56) { this->y56 = y56; }
  void SetThr(double thr) { this->thr = thr; }
  void SetHJM(double hjm) { this->hjm = hjm; }
  void SetLJM(double ljm) { this->ljm = ljm; }
  void SetWJB(double wjb) { this->wjb = wjb; }
  void SetNJB(double njb) { this->njb = njb; }

  void SetTAxis(Vec4 t_axis) { this->t_axis = t_axis; }

  // Member Functions ----------------------------------------------------------

  // Validation of Result Data
  bool Validate() {
    Vec4 psum = Vec4();

    std::vector<int> csum(100, 0);

    for (int i = 0; i < GetSize(); i++) {
      Parton p = GetParton(i);

      Vec4 pmom = p.GetMom();
      int pcol = p.GetCol();
      int pAntiCol = p.GetAntiCol();

      psum = psum + pmom;

      if (pcol > 0) {
        csum[pcol] += 1;
      }

      if (pAntiCol > 0) {
        csum[pAntiCol] -= 1;
      }
    }

    bool pcheck = (psum[0] < 1e-12 && psum[1] < 1e-12 && psum[2] < 1e-12 &&
                   psum[3] < 1e-12);
    if (!pcheck) {
      std::cout << psum << std::endl;
    }

    bool ccheck = true;
    for (int i = 0; i < 100; i++) {
      if (csum[i] != 0) {
        std::cout << "Colour " << i << " is not conserved." << std::endl;
        ccheck = false;
        break;
      }
    }

    return pcheck && ccheck;
  }

  void print_info() const {
    std::cout << "Event Information:\n";
    std::cout << "Dxs: " << GetDxs() << "\n";
    std::cout << "Number of Emissions: " << GetEmissions() << "\n";
    std::cout << "Shower T: " << GetShowerT() << "\n";
    std::cout << "Shower Y: " << GetShowerY() << "\n";
    std::cout << "Shower Z: " << GetShowerZ() << "\n";
    std::cout << "Shower C: " << GetShowerC() << "\n";
    std::cout << "Winner SF: " << GetWinSF() << "\n";
    std::cout << "Winner Dipole 1: " << GetWinDipole(0) << "\n";
    std::cout << "Winner Dipole 2: " << GetWinDipole(1) << "\n";
    std::cout << "Winner Params 1: " << GetWinParam(0) << "\n";
    std::cout << "Winner Params 2: " << GetWinParam(1) << "\n";
    std::cout << "Partons:\n";
    for (int i = 0; i < GetSize(); i++) {
      Parton parton = GetParton(i);
      std::cout << "  Parton " << i << ":\n";
      std::cout << "    Pid: " << parton.GetPid() << "\n";
      std::cout << "    Mom: " << parton.GetMom() << "\n";
      std::cout << "    Col: " << parton.GetCol() << "\n";
      std::cout << "    AntiCol: " << parton.GetAntiCol() << "\n";
    }
  }
};

#endif  // PARTON_H_
