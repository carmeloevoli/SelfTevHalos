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

}  // namespace CRWAVES