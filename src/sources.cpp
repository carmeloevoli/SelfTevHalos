// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow4 utils::pow_integer<4>

double Waves::compute_constant_CR_source_term_1D() {
  const double x_min = par.source_pmin / mks::electron_mass_c;
  const double x_cut = par.source_cutoff / mks::electron_mass_c;
  const double I = I_of_alpha(par.alpha, x_min, x_cut);
  const double R = par.tube_radius;
  double out = par.spin_down_luminosity * pow2(1. + par.age / par.source_tdecay);
  out /= 4.0 * pow2(M_PI) * mks::c_light * pow4(mks::electron_mass_c) * pow2(R) * I;
  return out;
}

double Waves::compute_constant_CR_source_term_3D() {
  const double x_min = par.source_pmin / mks::electron_mass_c;
  const double x_cut = par.source_cutoff / mks::electron_mass_c;
  const double I = I_of_alpha(par.alpha, x_min, x_cut);
  double out = par.spin_down_luminosity * pow2(1. + par.age / par.source_tdecay);
  out /= 4. * M_PI * mks::c_light * pow4(mks::electron_mass_c) * I;
  return out;
}

double Waves::compute_constant_CR_source_term() {
  return (par.do_3D) ? compute_constant_CR_source_term_3D() : compute_constant_CR_source_term_1D();
}

void Waves::build_CR_source_term() {
  const double q0 = compute_constant_CR_source_term();
  const double size = par.source_size;
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double F = pl_cutoff(p.at(ip), par.alpha, par.source_pmin, par.source_cutoff);
    for (size_t iz = 0; iz < z.size(); ++iz) {
      const double z_abs = fabs(z.at(iz));
      const double G = (par.do_3D) ? gaussian_3D(z_abs, size) : gaussian_1D(z_abs, size);
      Q_cr.get(ip, iz) = q0 * G * F;
    }
  }
  Q_cr.show_grid("Q_CR", 1.);
}

}  // namespace CRWAVES