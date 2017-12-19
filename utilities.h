#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_integration.h>
#include "units.h"

double beta(const double& p);

double larmor_radius(const double& p, const double& B);

double inverse_larmor_radius(const double& k, const double& B);

size_t find_nearest(const std::vector<double>& x, const double& x_0);

void print_vector(const std::vector<double>& x);

double spectrum(const double& p, const double& alpha, const double& p_min, const double& p_cutoff);

double source_profile_1D(const double& z, const double& size);

double source_profile_3D(const double& z, const double& size);

double I_of_alpha(double alpha, double x_min, double x_c);

double S_i(const double& T_i, const double& gamma_e);

#endif /* UTILITIES_H_ */
