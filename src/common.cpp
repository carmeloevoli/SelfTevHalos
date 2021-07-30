// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "common.h"

#include <gsl/gsl_integration.h>
#include <utilities/misc.h>

#include <cmath>

#include "units.h"

namespace CRWAVES {

double beta(const double& p) {
  using std::sqrt;
  const double gamma2 = 1. + utils::pow_integer<2>(p / mks::electron_mass_c);
  return sqrt(1. - 1. / gamma2);
}

double larmor_radius(const double& p, const double& B) { return p / mks::elementary_charge / B; }

double inverse_larmor_radius(const double& k, const double& B) { return mks::elementary_charge * B / k; }

double f(double x, void* params) {
  double alpha = *(double*)params;
  double x_cutoff = *((double*)params + 1);
  double f = pow(x, 3. - alpha) * exp(-utils::pow_integer<2>(x / x_cutoff));
  return f;
}

double I_of_alpha(const double& alpha, const double& x_min, const double& x_cutoff) {
  int LIMIT = 10000;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  double result, error;
  double params[2] = {alpha, x_cutoff};
  gsl_function F;
  F.function = &f;
  F.params = &params;
  gsl_integration_qag(&F, x_min, 2.0 * x_cutoff, 0, 1e-3, LIMIT, 3, w, &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

double source_evolution(const double& t_now, const double& t_decay) { return std::pow(1. + t_now / t_decay, -2.); }

double source_profile_1D(const double& z, const double& size) {
  using std::exp;
  using std::pow;
  return pow(2.0 * M_PI * utils::pow_integer<2>(size), -1. / 2.) * exp(-0.5 * utils::pow_integer<2>(z / size));
}

double source_profile_3D(const double& z, const double& size) {
  using std::exp;
  using std::pow;
  return pow(2.0 * M_PI * utils::pow_integer<2>(size), -3. / 2.) * exp(-0.5 * utils::pow_integer<2>(z / size));
}

double source_spectrum(const double& p, const double& alpha, const double& p_min, const double& p_cutoff) {
  using std::exp;
  using std::pow;
  return (p > p_min) ? pow(p / mks::electron_mass_c, -alpha) * exp(-utils::pow_integer<2>(p / p_cutoff)) : 0.;
}

}  // namespace CRWAVES