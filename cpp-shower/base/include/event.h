#ifndef EVENT_H_
#define EVENT_H_

#include "parton.h"

// Event Class
// Built to contain the partons and the dxs as one accesible object
// In future, can use to store thrust, log10y23 to parallelise those
class Event {
 private:
  // Temporary Solution - Allows a limited number of partons
  // Better Solution would be to use a dynamic array, but not GPU friendly
  Parton partons[maxPartons];

  // ME Params -----------------------------------------------------------------

  double dxs = 0.;  // Differential Cross Section
  int nHard = 0;    // Number of Hard Partons
  // int nInitial = 0;    // Number of Initial Partons (Prep for ISR)
  // int nNonParton = 0;  // Number of Non-Parton Partons (Prep for ISR)

  // Shower Params -------------------------------------------------------------

  int nEmission = 0;    // Number of Emissions
  double showerT = 0.;  // Evolution and Splitting Variables
  double showerZ = 0.;
  double showerY = 0.;
  int showerC = 0;  // Colour Counter

  // Selecting Winner Emission - Defaults Values which represent no winner
  int winSF = 16;
  int winDipole[2] = {-1, -1};
  double winParams[2] = {0., 0.};

  bool endShower = false;  // Shower End Flag - used if T < 1 GeV

  // Analysis Variables --------------------------------------------------------

  // Event Validity - Momentum and Colour Conservation
  bool validity = true;

  // Jet Rates using the Durham Algorithm
  double y23 = -50., y34 = -50., y45 = -50., y56 = -50.;

  // Event Shape Variables - Thrust, Jet Masses and Broadenings
  double thr = -50., hjm = -50., ljm = -50., wjb = -50., njb = -50.;
  Vec4 t_axis = Vec4();

  // Dalitz Plot
  double dalitz[2] = {-50., -50.};

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
  int GetPartonSize() const { return (nHard + nEmission) - 2; }  // -2: e+, e-

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

  double GetDalitz(int i) const { return dalitz[i]; }

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

  void SetDalitz(double x1, double x2) {
    dalitz[0] = x1;
    dalitz[1] = x2;
  }

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

#endif  // EVENT_H_