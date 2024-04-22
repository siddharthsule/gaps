#include "eventshapes.h"

// Event Shapes

void bubbleSort(Vec4* moms, int n) {
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - i - 1; j++) {
      if (moms[j].P() < moms[j + 1].P()) {
        Vec4 temp = moms[j];
        moms[j] = moms[j + 1];
        moms[j + 1] = temp;
      }
    }
  }
}

void CalculateThrust(Event& ev) {
  if (!ev.GetValidity() || ev.GetPartonSize() < 3) {
    return;
  }

  Vec4 moms[maxPartons];
  for (int i = 2; i < ev.GetSize(); ++i) {
    moms[i - 2] = ev.GetParton(i).GetMom();
  }

  bubbleSort(moms, maxPartons);

  double momsum = 0.;
  for (int i = 0; i < ev.GetPartonSize(); ++i) {
    momsum += moms[i].P();
  }

  double thr = 0.;
  Vec4 t_axis = Vec4();

  for (int k = 1; k < ev.GetPartonSize(); ++k) {
    for (int j = 0; j < k; ++j) {
      Vec4 tmp_axis = moms[j].Cross(moms[k]);
      Vec4 p_thrust = Vec4();
      Vec4 p_combin[4];

      for (int i = 0; i < ev.GetPartonSize(); ++i) {
        if (i != j && i != k) {
          if (moms[i].Dot(tmp_axis) >= 0) {
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
        double temp = p_combin[i].P();
        if (temp > thr) {
          thr = temp;
          t_axis = p_combin[i];
        }
      }
    }
  }

  thr /= momsum;
  thr = 1. - thr;

  t_axis = t_axis / (t_axis).P();
  if (t_axis[3] < 0) {
    t_axis = t_axis * -1.;
  }

  if (thr < 1e-12) {
    thr = -5.;
  }

  ev.SetThr(thr);
  ev.SetTAxis(t_axis);
}

void CalculateJetMBr(Event& ev) {
  if (!ev.GetValidity() || ev.GetPartonSize() < 3) {
    return;
  }

  Vec4 moms[maxPartons];
  for (int i = 2; i < ev.GetSize(); ++i) {
    moms[i - 2] = ev.GetParton(i).GetMom();
  }

  double momsum = 0.;
  for (int i = 0; i < ev.GetSize(); ++i) {
    momsum += moms[i].P();
  }

  Vec4 p_with, p_against;
  int n_with = 0, n_against = 0;
  double e_vis = 0., broad_with = 0., broad_against = 0.,
         broad_denominator = 0.;

  for (int i = 0; i < ev.GetPartonSize(); ++i) {
    double mo_para = moms[i].Dot(ev.GetTAxis());
    double mo_perp = (moms[i] - (ev.GetTAxis() * mo_para)).P();
    double enrg = moms[i].P();

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

  double mass2_with = fabs(p_with.M2() / e2_vis);
  double mass2_against = fabs(p_against.M2() / e2_vis);

  double mass_with = sqrt(mass2_with);
  double mass_against = sqrt(mass2_against);

  broad_with /= broad_denominator;
  broad_against /= broad_denominator;

  double mH = fmax(mass_with, mass_against);
  double mL = fmin(mass_with, mass_against);

  double bW = fmax(broad_with, broad_against);
  double bN = fmin(broad_with, broad_against);

  if (n_with == 1 || n_against == 1) {
    ev.SetHJM(mH);
    ev.SetWJB(bW);
  } else {
    ev.SetHJM(mH);
    ev.SetLJM(mL);
    ev.SetWJB(bW);
    ev.SetNJB(bN);
  }
}