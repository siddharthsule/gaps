#include "particle.cuh"

// Boolean - If two partons are in a Colour Connected Pair
bool Particle::IsColorConnected(Particle p) {
  return (col > 0 && col == p.anticol) || (anticol > 0 && anticol == p.col);
}


// Validation of Parton Shower Result sData
bool IsEventCheckValid(const Event& ev) {
  Vec4 psum = Vec4();

  std::vector<int> csum(100, 0);

  for (int i = 0; i < ev.GetSize(); i++) {
    
    Particle p = ev.GetParton(i);
    
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

  bool pcheck = (psum[0] < 1e-12 && psum[1] < 1e-12 && psum[2] < 1e-12 && psum[3] < 1e-12);
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

