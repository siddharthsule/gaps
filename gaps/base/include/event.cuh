#ifndef EVENT_CUH_
#define EVENT_CUH_

#include "parton.cuh"

/**
 * The Event Class
 * ---------------
 *
 * This is the backbone of the program. An event contains the differential cross
 * section calculated from the ME, the hard partons and the showered partons.
 * It also contains the shower parameters which are constantly updated during
 * the showering process. The event also contains the analysis variables, which
 * are calculated after the showering process is complete.
 */

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

  // End Shower Flag
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
  __device__ Event() {}

  // Getters -------------------------------------------------------------------

  // Access Partons in the Event
  __device__ Parton GetParton(int i) const { return partons[i]; }
  __device__ int GetSize() const { return nHard + nEmission; }
  __device__ int GetHard() const { return nHard; }
  __device__ int GetEmissions() const { return nEmission; }
  __device__ int GetPartonSize() const {
    // -2: e+, e-
    return (nHard + nEmission) - 2;
  }

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

  __device__ double GetDalitz(int i) const { return dalitz[i]; }

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

  __device__ void SetDalitz(double x1, double x2) {
    dalitz[0] = x1;
    dalitz[1] = x2;
  }

  // Member Functions ----------------------------------------------------------

  // Validate the Event - Check Momentum and Colour Conservation
  __device__ bool Validate() {
    Vec4 psum = Vec4();

    // N Colours = N Partons - 1
    int csum[maxPartons - 1] = {0};

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

    bool pcheck = (abs(psum[0]) < 1e-12 && abs(psum[1]) < 1e-12 &&
                   abs(psum[2]) < 1e-12 && abs(psum[3]) < 1e-12);

    /* // No need to print for GPU, it counts number of invalid events
    if (!pcheck) {
      printf("%f %f %f %f\n", psum[0], psum[1], psum[2], psum[3]);
    }
    */

    bool ccheck = true;
    for (int i = 0; i < maxPartons - 1; i++) {
      if (csum[i] != 0) {
        // printf("Colour %d is not conserved.\n", i);
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

#endif  // EVENT_CUH_