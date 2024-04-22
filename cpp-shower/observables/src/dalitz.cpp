#include "dalitz.h"

// Dalitz Plot

void CalculateDalitz(Event& ev) {
  if (!ev.GetValidity() || ev.GetPartonSize() != 3) {
    return;
  }

  // Obtain Energy from incoming partons
  double E = abs(ev.GetParton(0).GetMom()[0] + ev.GetParton(1).GetMom()[0]);

  // By default, element 2 is quark and 3 is antiquark
  // i.e. emission will be element 4
  Vec4 p1 = ev.GetParton(2).GetMom();
  Vec4 p2 = ev.GetParton(3).GetMom();

  // Calculate x1 and x2
  double x1 = 2 * p1.P() / E;
  double x2 = 2 * p2.P() / E;

  ev.SetDalitz(x1, x2);
}