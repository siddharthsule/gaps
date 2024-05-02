#include "eventshapes.h"

// event shapes

void bubble_sort(vec4* moms, int n) {
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

void calculate_thrust(event& ev) {
  if (!ev.get_validity() || ev.get_parton_size() < 3) {
    return;
  }

  vec4 moms[max_partons];
  for (int i = 2; i < ev.get_size(); ++i) {
    moms[i - 2] = ev.get_parton(i).get_mom();
  }

  bubble_sort(moms, max_partons);

  double momsum = 0.;
  for (int i = 0; i < ev.get_parton_size(); ++i) {
    momsum += moms[i].p();
  }

  double thr = 0.;
  vec4 t_axis = vec4();

  for (int k = 1; k < ev.get_parton_size(); ++k) {
    for (int j = 0; j < k; ++j) {
      vec4 tmp_axis = moms[j].cross(moms[k]);
      vec4 p_thrust = vec4();
      vec4 p_combin[4];

      for (int i = 0; i < ev.get_parton_size(); ++i) {
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

  if (thr < 1e-12) {
    thr = -5.;
  }

  ev.set_thr(thr);
  ev.set_t_axis(t_axis);
}

void calculate_jet_m_br(event& ev) {
  if (!ev.get_validity() || ev.get_parton_size() < 3) {
    return;
  }

  vec4 moms[max_partons];
  for (int i = 2; i < ev.get_size(); ++i) {
    moms[i - 2] = ev.get_parton(i).get_mom();
  }

  double momsum = 0.;
  for (int i = 0; i < ev.get_size(); ++i) {
    momsum += moms[i].p();
  }

  vec4 p_with, p_against;
  int n_with = 0, n_against = 0;
  double e_vis = 0., broad_with = 0., broad_against = 0.,
         broad_denominator = 0.;

  for (int i = 0; i < ev.get_parton_size(); ++i) {
    double mo_para = moms[i].dot(ev.get_t_axis());
    double mo_perp = (moms[i] - (ev.get_t_axis() * mo_para)).p();
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

  if (n_with == 1 || n_against == 1) {
    ev.set_hjm(m_h);
    ev.set_wjb(b_w);
  } else {
    ev.set_hjm(m_h);
    ev.set_ljm(m_l);
    ev.set_wjb(b_w);
    ev.set_njb(b_n);
  }
}