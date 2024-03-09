#ifndef PARTON_CUH_
#define PARTON_CUH_

// Partons have Vec4 Momentum, Vec4 #includes Base
#include "vec4.cuh"

/**
 * The Parton and Event Classes
 * ----------------------------
 *
 * This file is the backbone of the entire program. It contains the Parton and
 * Event classes, which are used to store the partons and the event.
 *
 * An event contains the differential cross section calculated from the ME, the
 * hard partons and the showered partons. It also contains the shower parameters
 * which are constantly updated during the showering process. The event also
 * contains the analysis variables, which are calculated after the showering
 * process is complete.
 */

class Parton {
 public:
  // Constructor
  __device__ Parton(int pid = 0, Vec4 momentum = Vec4(), int col = 0,
                    int anticol = 0)
      : pid(pid), mom(momentum), col(col), anticol(anticol) {}

  // Getters and Setters
  __device__ int GetPid() const { return pid; }
  __device__ Vec4 GetMom() const { return mom; }
  __device__ int GetCol() const { return col; }
  __device__ int GetAntiCol() const { return anticol; }

  __device__ void SetPid(int pid) { this->pid = pid; }
  __device__ void SetMom(Vec4 mom) { this->mom = mom; }
  __device__ void SetCol(int col) { this->col = col; }
  __device__ void SetAntiCol(int anticol) { this->anticol = anticol; }

  // If two partons are in a Colour Connected Dipole
  __device__ bool IsColorConnected(Parton p) {
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
class Event {
 private:
  // Temporary Solution - Allows a limited number of partons
  // Better Solution would be to use a dynamic array, but not GPU friendly
  Parton partons[maxPartons];

  // ME Params -----------------------------------------------------------------

  double dxs; // Differential Cross Section
  int nHard; // Number of Hard Partons

  // Shower Params -------------------------------------------------------------

  int nEmission; // Number of Emissions
  double showerT, showerZ, showerY; // Evolution and Splitting Variables
  int showerC; // Colour Counter

  // Selecting Winner Emission - Defaults Values which represent no winner
  int winSF = 16;
  int winDipole[2] = {-1, -1};
  double winParams[2] = {0.0, 0.0};

  bool endShower = false; // Shower End Flag - used if T < 1 GeV

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
  __device__ Event() {}

  // Getters -------------------------------------------------------------------

  // Access Partons in the Event
  __device__ Parton GetParton(int i) const { return partons[i]; }
  __device__ int GetSize() const { return nHard + nEmission; }
  __device__ int GetHard() const { return nHard; }
  __device__ int GetEmissions() const { return nEmission; }
  __device__ int GetPartonSize() const { return nHard + nEmission - 2; }

  // Get Differential Cross Section
  __device__ double GetDxs() const { return dxs; }

  // Get Shower Params
  __device__ double GetShowerT() const { return showerT; }
  __device__ double GetShowerZ() const { return showerZ; }
  __device__ double GetShowerY() const { return showerY; }
  __device__ int GetShowerC() const { return showerC; }

  __device__ int GetWinSF() const { return winSF; }
  __device__ int GetWinDipole(int i) const { return winDipole[i]; }
  __device__ double GetWinParam(int i) const { return winParams[i]; }

  __device__ bool GetEndShower() const { return endShower; }

  // Analysis Getters
  __device__ bool GetValidity() const { return validity; }

  __device__ double GetY23() const { return y23; }
  __device__ double GetY34() const { return y34; }
  __device__ double GetY45() const { return y45; }
  __device__ double GetY56() const { return y56; }
  __device__ double GetThr() const { return thr; }
  __device__ double GetHJM() const { return hjm; }
  __device__ double GetLJM() const { return ljm; }
  __device__ double GetWJB() const { return wjb; }
  __device__ double GetNJB() const { return njb; }

  __device__ Vec4 GetTAxis() const { return t_axis; }

  // Setters -------------------------------------------------------------------
  
  // Add / Replace Parton
  __device__ void SetParton(int i, Parton parton) { partons[i] = parton; }
  
  // Set Parton Data
  __device__ void SetPartonPid(int i, int pid) { partons[i].SetPid(pid); }
  __device__ void SetPartonMom(int i, Vec4 mom) { partons[i].SetMom(mom); }
  __device__ void SetPartonCol(int i, int col) { partons[i].SetCol(col); }
  __device__ void SetPartonAntiCol(int i, int anticol) {
    partons[i].SetAntiCol(anticol);
  }

  // Set Differential Cross Section and nHard
  __device__ void SetDxs(double dxs) { this->dxs = dxs; }
  __device__ void SetHard(int nHard) { this->nHard = nHard; }

  // Adjust and Increment Number of Emissions
  __device__ void SetEmissions(int nEmission) { this->nEmission = nEmission; }
  __device__ void IncrementEmissions() { nEmission++; }

  // Set Shower Params
  __device__ void SetShowerT(double showerT) { this->showerT = showerT; }
  __device__ void SetShowerZ(double showerZ) { this->showerZ = showerZ; }
  __device__ void SetShowerY(double showerY) { this->showerY = showerY; }

  __device__ void SetShowerC(int showerC) { this->showerC = showerC; }
  __device__ void IncrementShowerC() { showerC++; }

  __device__ void SetWinSF(int winSF) { this->winSF = winSF; }
  __device__ void SetWinDipole(int i, int winParton) {
    this->winDipole[i] = winParton;
  }
  __device__ void SetWinParam(int i, double winParam) {
    this->winParams[i] = winParam;
  }

  __device__ void SetEndShower(bool endShower) { this->endShower = endShower; }

  // Set Analysis Variables
  __device__ void SetValidity(bool validity) { this->validity = validity; }

  __device__ void SetY23(double y23) { this->y23 = y23; }
  __device__ void SetY34(double y34) { this->y34 = y34; }
  __device__ void SetY45(double y45) { this->y45 = y45; }
  __device__ void SetY56(double y56) { this->y56 = y56; }
  __device__ void SetThr(double thr) { this->thr = thr; }
  __device__ void SetHJM(double hjm) { this->hjm = hjm; }
  __device__ void SetLJM(double ljm) { this->ljm = ljm; }
  __device__ void SetWJB(double wjb) { this->wjb = wjb; }
  __device__ void SetNJB(double njb) { this->njb = njb; }

  __device__ void SetTAxis(Vec4 t_axis) { this->t_axis = t_axis; }

  // Member Functions ----------------------------------------------------------

  // Validate the Event - Check Momentum and Colour Conservation
  __device__ bool Validate() {
    Vec4 psum = Vec4();

    int csum[100] = {0};

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
      printf("%f %f %f %f\n", psum[0], psum[1], psum[2], psum[3]);
    }

    bool ccheck = true;
    for (int i = 0; i < maxPartons - 1; i++) {
      if (csum[i] != 0) {
        printf("Colour %d is not conserved.\n", i);
        ccheck = false;
        break;
      }
    }

    return pcheck && ccheck;
  }

  __device__ void print_info() const {
    printf("Event Information:\n");
    printf("Dxs: %f\n", GetDxs());
    printf("Number of Emissions: %d\n", GetEmissions());
    printf("Shower Parameters:\n");
    printf("  T: %f\n", GetShowerT());
    printf("  Y: %f\n", GetShowerY());
    printf("  Z: %f\n", GetShowerZ());
    printf("  C: %d\n", GetShowerC());
    printf("Shower Winner:\n");
    printf("  Kernel Number: %d\n", GetWinSF());
    printf("  Partons: [%d, %d]\n", GetWinDipole(0), GetWinDipole(1));
    printf("  Params: [%f, %f]\n", GetWinParam(0), GetWinParam(1));
    printf("Partons:\n");
    for (int i = 0; i < GetSize(); i++) {
      Parton parton = GetParton(i);
      printf("  Parton %d:\n", i);
      printf("    Pid: %d\n", parton.GetPid());
      printf("    Mom: %f\n", parton.GetMom().P());
      printf("    Col: %d\n", parton.GetCol());
      printf("    AntiCol: %d\n", parton.GetAntiCol());
    }
  }
};

#endif  // PARTON_CUH_