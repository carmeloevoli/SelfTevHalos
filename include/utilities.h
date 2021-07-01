// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <gsl/gsl_integration.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "units.h"

namespace CRWAVES {

double beta(const double& p);

double larmor_radius(const double& p, const double& B);

double inverse_larmor_radius(const double& k, const double& B);

size_t find_nearest(const std::vector<double>& x, const double& x_0);

void print_vector(const std::vector<double>& x);

double source_evolution(const double& t_now, const double& t_decay);

}  // namespace CRWAVES

#endif /* UTILITIES_H_ */
