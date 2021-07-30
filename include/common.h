// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <utilities/misc.h>

#include <cmath>

#include "units.h"

namespace CRWAVES {

inline double larmor_radius(double p, double B) { return p / mks::elementary_charge / B; }
inline double inverse_larmor_radius(double k, double B) { return mks::elementary_charge * B / k; }
inline double source_evolution(double t, double t_decay) { return 1. / utils::pow_integer<2>(1. + t / t_decay); }
inline double gamma2(double p) { return 1. + utils::pow_integer<2>(p / mks::electron_mass_c); }
inline double beta(double p) { return std::sqrt(1. - 1. / gamma2(p)); }
inline double gaussian_1D(double z, double size) {
  return std::pow(2.0 * M_PI * utils::pow_integer<2>(size), -0.5) * std::exp(-0.5 * utils::pow_integer<2>(z / size));
}
inline double gaussian_3D(double z, double size) {
  return std::pow(2.0 * M_PI * utils::pow_integer<2>(size), -1.5) * std::exp(-0.5 * utils::pow_integer<2>(z / size));
}
inline double pl_cutoff(const double& p, const double& alpha, const double& p_min, const double& p_cutoff) {
  return (p > p_min) ? std::pow(p / mks::electron_mass_c, -alpha) * std::exp(-utils::pow_integer<2>(p / p_cutoff)) : 0.;
}

double I_of_alpha(const double& alpha, const double& x_min, const double& x_cutoff);

}  // namespace CRWAVES

#endif /* UTILITIES_H_ */
