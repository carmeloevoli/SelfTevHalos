// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "common.h"

#include <gsl/gsl_integration.h>

namespace CRWAVES {

double f(double x, void* params) {
  double alpha = *(double*)params;
  double x_cutoff = *((double*)params + 1);
  double f = std::pow(x, 3. - alpha) * std::exp(-(x / x_cutoff));
  return f;
}

double pl_integral(const double& alpha, const double& x_min, const double& x_cutoff) {
  int LIMIT = 10000;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  double result, error;
  double params[2] = {alpha, x_cutoff};
  gsl_function F;
  F.function = &f;
  F.params = &params;
  gsl_integration_qag(&F, x_min, 2.0 * x_cutoff, 0, 1e-5, LIMIT, 3, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double KN_factor(const double& T, const double& gamma_e) {
  double value = 45. * pow2(cgs::me_c2) / 64. / pow2(M_PI) / pow2(cgs::k_B * T);
  return value / (value + pow2(gamma_e));
}

double momentum_loss_rate(double gamma_e, double B) {
  double energy_density = 0.26 * cgs::eV / cgs::cm3 * KN_factor(2.7 * cgs::K, gamma_e);  // CMB
  energy_density += 0.30 * cgs::eV / cgs::cm3 * KN_factor(20 * cgs::K, gamma_e);         // IR
  energy_density += 0.30 * cgs::eV / cgs::cm3 * KN_factor(5000 * cgs::K, gamma_e);       // star
  energy_density += 0.10 * cgs::eV / cgs::cm3 * KN_factor(20000 * cgs::K, gamma_e);      // UV
  energy_density += magnetic_energy_density(B);
  return -4. / 3. * cgs::sigma_th * energy_density * pow2(gamma_e);
}

double D2W(double D, double p, double B) {
  const double r_L = larmor_radius(p, B);
  return cgs::c_light * utils::pow_integer<2>(r_L) / 3. / D;
}

double W2D(double W, double p, double B) {
  const double r_L = larmor_radius(p, B);
  return cgs::c_light * utils::pow_integer<2>(r_L) / 3. / W;
}

double source_term_norm_1D(double p_min, double p_cutoff, double slope, double L_0, double R) {
  const double x_min = p_min / cgs::me_c;
  const double x_cut = p_cutoff / cgs::me_c;
  const double I = pl_integral(slope, x_min, x_cut);
  double value = L_0;
  value /= 4.0 * cgs::c_light * utils::pow_integer<4>(cgs::me_c) * utils::pow_integer<2>(M_PI * R) * I;
  return value;
}

double source_term_norm_3D(double p_min, double p_cutoff, double slope, double L_0, double R) {
  //   const double x_min = par.source_pmin / cgs::me_c;
  //   const double x_cut = par.source_cutoff / cgs::me_c;
  //   const double I = pl_integral(par.source_slope, x_min, x_cut);
  //   double value = initial_luminosity(par.source_luminosity_today, par.source_age, par.source_tdecay);
  //   value /= 4. * M_PI * cgs::c_light * pow4(cgs::me_c) * I;
  //   return value;
  return 0;
}

}  // namespace CRWAVES