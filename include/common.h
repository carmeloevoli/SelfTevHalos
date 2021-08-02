// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cmath>

#include "utilities/misc.h"
#include "utilities/units.h"

namespace CRWAVES {

inline double larmor_radius(double p, double B) { return p * cgs::c_light / cgs::elementary_charge / B; }
inline double resonant_momentum(double k, double B) { return cgs::elementary_charge * B / cgs::c_light / k; }
inline double alfvenSpeed(double ion_density, double B) {
  return B / std::sqrt(4 * M_PI * cgs::proton_mass * ion_density);
}
inline double magnetic_energy_density(double B) { return B * B / 8. / M_PI; }
inline double source_evolution(double t, double t_decay) { return 1. / utils::pow_integer<2>(1. + t / t_decay); }
inline double initial_luminosity(double L_today, double age, double t_decay) {
  return L_today * utils::pow_integer<2>(1. + age / t_decay);
}
inline double gamma2(double p) { return 1. + utils::pow_integer<2>(p / cgs::me_c); }
inline double beta(double p) { return std::sqrt(1. - 1. / gamma2(p)); }
inline double gaussian_1D(double z, double size) {
  return std::pow(2.0 * M_PI * utils::pow_integer<2>(size), -0.5) * std::exp(-0.5 * utils::pow_integer<2>(z / size));
}
inline double gaussian_3D(double z, double size) {
  return std::pow(2.0 * M_PI * utils::pow_integer<2>(size), -1.5) * std::exp(-0.5 * utils::pow_integer<2>(z / size));
}
inline double pl_with_cutoff(const double& p, const double& alpha, const double& p_min, const double& p_cutoff) {
  return (p > p_min) ? std::pow(p / cgs::me_c, -alpha) * std::exp(-(p / p_cutoff)) : 0.;
}

double pl_integral(const double& alpha, const double& x_min, const double& x_cutoff);

}  // namespace CRWAVES

#endif /* UTILITIES_H_ */
