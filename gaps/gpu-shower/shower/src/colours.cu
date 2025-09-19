#include "shower.cuh"

__device__ void shower::make_colours(int current_col, int sf, int flavs[3],
                                     int colij[2], int colk[2], int* coli,
                                     int* colj, double r) const {
  /**
   * @brief Generate the final colour state for the given splitting
   *
   * @param current_col: Current number of colours
   * @param sf: Splitting Function Code
   * @param flavs: Flavour array
   * @param colij: Colour of the emitter
   * @param colk: Colour of the spectator
   * @param coli: Colour of the emitter after the splitting
   * @param colj: Colour of the emitted parton
   * @param r: Random number between 0 and 1
   */

  // New Colour = Current Number of Colours + 1
  int new_col = current_col + 1;

  // Final State Emission
  if (is_fsr(sf)) {
    // q -> qg
    if (flavs[0] != 21) {
      // q emitter
      if (flavs[0] > 0) {
        coli[0] = new_col;
        coli[1] = 0;
        colj[0] = colij[0];
        colj[1] = new_col;
      }
      // qbar -> qbar g
      else {
        coli[0] = 0;
        coli[1] = new_col;
        colj[0] = new_col;
        colj[1] = colij[1];
      }
    }

    else {
      // g -> gg
      if (flavs[1] == 21) {
        if (colij[0] == colk[1]) {
          // Criss-cross for Colour Single Gluon Dipole
          if (colij[1] == colk[0] && r > 0.5) {
            coli[0] = colij[0];
            coli[1] = new_col;
            colj[0] = new_col;
            colj[1] = colij[1];
          } else {
            coli[0] = new_col;
            coli[1] = colij[1];
            colj[0] = colij[0];
            colj[1] = new_col;
          }
        } else {
          coli[0] = colij[0];
          coli[1] = new_col;
          colj[0] = new_col;
          colj[1] = colij[1];
        }
      }
      // g -> qq
      else {
        // g -> q qbar
        if (flavs[1] > 0) {
          coli[0] = colij[0];
          coli[1] = 0;
          colj[0] = 0;
          colj[1] = colij[1];
        }
        // g -> qbar q
        else {
          coli[0] = 0;
          coli[1] = colij[1];
          colj[0] = colij[0];
          colj[1] = 0;
        }
      }
    }
  }

  // Initial State Emission - Incoming Parton has Col <-> AntiCol
  else {
    // q emitter
    if (flavs[0] != 21) {
      // q -> qg
      if (flavs[1] != 21) {
        if (flavs[0] > 0) {
          // Incoming, so quark updates anticolour
          coli[0] = 0;
          coli[1] = new_col;
          // Emitted gluon gets new colour and quark's anticolour
          colj[0] = new_col;
          colj[1] = colij[1];
        } else {
          // Incoming, so antiquark updates colour
          coli[0] = new_col;
          coli[1] = 0;
          // Emitted gluon gets antiquark's colour and new anticolour
          colj[0] = colij[0];
          colj[1] = new_col;
        }
      }
      // q -> gq (backwards is g -> qq)
      else {
        if (flavs[0] > 0) {
          // Incoming: gluon gets new colour and quark's anticolour
          coli[0] = new_col;
          coli[1] = colij[1];
          // Emitted (anti)quark gets new anticolour
          colj[0] = 0;
          colj[1] = new_col;
        } else {
          // Incoming: gluon gets antiquark's colour and new anticolour
          coli[0] = colij[0];
          coli[1] = new_col;
          // Emitted quark gets new colour
          colj[0] = new_col;
          colj[1] = 0;
        }
      }
    }

    else {
      // g -> gg - criss-cross is not colour supressed
      if (flavs[1] == 21) {
        if (colij[0] == colk[1]) {
          coli[0] = new_col;
          coli[1] = colij[1];
          colj[0] = colij[0];
          colj[1] = new_col;
        } else {
          coli[0] = colij[0];
          coli[1] = new_col;
          colj[0] = new_col;
          colj[1] = colij[1];
        }
      }

      // g -> qq (backwards is q -> gq)
      else {
        // g -> qbar q (backwards is q -> g q)
        if (flavs[1] > 0) {
          // Incoming: quark gets current anticolour
          coli[0] = 0;
          coli[1] = colij[1];
          // Emitted quark gets current colour
          colj[0] = colij[0];
          colj[1] = 0;
        }
        // g -> q qbar (backwards is qbar -> g qbar)
        else {
          // Incoming: antiquark gets current colour
          coli[0] = colij[0];
          coli[1] = 0;
          // Emitted antiquark gets current anticolour
          colj[0] = 0;
          colj[1] = colij[1];
        }
      }
    }
  }
}