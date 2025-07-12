#include "eventshapes.cuh"

// event shapes

__device__ void bubble_sort(vec4* moms, int n) {
  /**
   * @brief Bubble Sort the Momenta in Descending Order of p()
   *
   * @param moms The array of momenta
   * @param n The number of momenta
   */

  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - i - 1; j++) {
      if (moms[j].p() < moms[j + 1].p()) {
        vec4 temp = moms[j];
        moms[j] = moms[j + 1];
        moms[j + 1] = temp;
      }
    }
  }
}

__global__ void calculate_ev_shapes(const event* events, double* results,
                                    int n) {
  /**
   * @brief Calculate the Event Shapes for LEP
   *
   * @param ev The event object
   * @param results The array to store the results
   * @param n The number of events
   */
  // ---------------------------------------------
  // Kernel Preamble
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx >= n) return;
  // ---------------------------------------------
  // Observables Preamble
  if (!events[idx].get_validity()) return;
  const event& ev = events[idx];
  // ---------------------------------------------

  // Event Shapes Limited to More than 2 particles
  if (ev.get_size() - 2 < 3) {
    return;
  }

  // Prepare the Momenta
  vec4 moms[max_particles];
  for (int i = 2; i < ev.get_size(); ++i) {
    moms[i - 2] = ev.get_particle(i).get_mom();
  }

  // Sort the Momenta
  bubble_sort(moms, max_particles);

  // Calculate Total Momentum
  double momsum = 0.;
  for (int i = 0; i < ev.get_size() - 2; ++i) {
    momsum += moms[i].p();
  }

  // --------------------------------------------------------------
  // Thrust Calculation, from TASSO

  double thr = 0.;
  vec4 t_axis = vec4();

  for (int k = 1; k < ev.get_size() - 2; ++k) {
    for (int j = 0; j < k; ++j) {
      vec4 tmp_axis = moms[j].cross(moms[k]);
      vec4 p_thrust = vec4();
      vec4 p_combin[4];

      for (int i = 0; i < ev.get_size() - 2; ++i) {
        if (i != j && i != k) {
          if (moms[i].dot(tmp_axis) >= 0) {
            p_thrust = p_thrust + moms[i];
          } else {
            p_thrust = p_thrust - moms[i];
          }
        }
      }

      p_combin[0] = (p_thrust + moms[j] + moms[k]);
      p_combin[1] = (p_thrust + moms[j] - moms[k]);
      p_combin[2] = (p_thrust - moms[j] + moms[k]);
      p_combin[3] = (p_thrust - moms[j] - moms[k]);

      for (int i = 0; i < 4; ++i) {
        double temp = p_combin[i].p();
        if (temp > thr) {
          thr = temp;
          t_axis = p_combin[i];
        }
      }
    }
  }

  thr /= momsum;
  thr = 1. - thr;

  t_axis = t_axis / (t_axis).p();
  if (t_axis[3] < 0) {
    t_axis = t_axis * -1.;
  }

  // DISABLE FOR LOG
  if (thr < 1e-12) {
    thr = -5.;
  }

  // --------------------------------------------------------------
  // Jet Mass and Broadening Calculation

  vec4 p_with, p_against;
  int n_with = 0, n_against = 0;
  double e_vis = 0., broad_with = 0., broad_against = 0.,
         broad_denominator = 0.;

  for (int i = 0; i < ev.get_size() - 2; ++i) {
    double mo_para = moms[i].dot(t_axis);
    double mo_perp = (moms[i] - (t_axis * mo_para)).p();
    double enrg = moms[i].p();

    e_vis += enrg;
    broad_denominator += 2. * enrg;

    if (mo_para > 0.) {
      p_with = p_with + moms[i];
      broad_with += mo_perp;
      n_with++;
    } else if (mo_para < 0.) {
      p_against = p_against + moms[i];
      broad_against += mo_perp;
      n_against++;
    } else {
      p_with = p_with + (moms[i] * 0.5);
      p_against = p_against + (moms[i] * 0.5);
      broad_with += 0.5 * mo_perp;
      broad_against += 0.5 * mo_perp;
      n_with++;
      n_against++;
    }
  }

  double e2_vis = e_vis * e_vis;

  double mass2_with = fabs(p_with.m2() / e2_vis);
  double mass2_against = fabs(p_against.m2() / e2_vis);

  double mass_with = sqrt(mass2_with);
  double mass_against = sqrt(mass2_against);

  broad_with /= broad_denominator;
  broad_against /= broad_denominator;

  double m_h = fmax(mass_with, mass_against);
  double m_l = fmin(mass_with, mass_against);

  double b_w = fmax(broad_with, broad_against);
  double b_n = fmin(broad_with, broad_against);

  // store the results (y23, y34, y45, y56, tvalue, tzoomd, hjm, ljm, wjb, njb)
  results[21 * idx + 4] = thr;
  results[21 * idx + 5] = thr;
  results[21 * idx + 6] = m_h;
  results[21 * idx + 7] = (n_with == 1 || n_against == 1) ? -50. : m_l;
  results[21 * idx + 8] = b_w;
  results[21 * idx + 9] = (n_with == 1 || n_against == 1) ? -50. : b_n;

  // Log Observables - for NLL Tests
  // results[21 * idx + 4] = log(thr);
  // results[21 * idx + 8] = log(b_w);
}