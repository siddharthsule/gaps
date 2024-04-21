#include "shower.h"

void Shower::MakeColours(Event& ev, int* coli, int* colj, const int flavs[3],
                         const int colij[2], const int colk[2],
                         const double r) {
  // Increase variable ev.GetShowerC() by 1
  ev.IncrementShowerC();

  if (flavs[0] != 21) {
    if (flavs[0] > 0) {
      coli[0] = ev.GetShowerC();
      coli[1] = 0;
      colj[0] = colij[0];
      colj[1] = ev.GetShowerC();
    } else {
      coli[0] = 0;
      coli[1] = ev.GetShowerC();
      colj[0] = ev.GetShowerC();
      colj[1] = colij[1];
    }
  } else {
    if (flavs[1] == 21) {
      if (colij[0] == colk[1]) {
        if (colij[1] == colk[0] && r > 0.5) {
          coli[0] = colij[0];
          coli[1] = ev.GetShowerC();
          colj[0] = ev.GetShowerC();
          colj[1] = colij[1];
        } else {
          coli[0] = ev.GetShowerC();
          coli[1] = colij[1];
          colj[0] = colij[0];
          colj[1] = ev.GetShowerC();
        }
      } else {
        coli[0] = colij[0];
        coli[1] = ev.GetShowerC();
        colj[0] = ev.GetShowerC();
        colj[1] = colij[1];
      }
    } else {
      if (flavs[1] > 0) {
        coli[0] = colij[0];
        coli[1] = 0;
        colj[0] = 0;
        colj[1] = colij[1];
      } else {
        coli[0] = 0;
        coli[1] = colij[1];
        colj[0] = colij[0];
        colj[1] = 0;
      }
    }
  }
}