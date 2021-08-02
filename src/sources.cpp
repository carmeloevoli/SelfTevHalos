// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow3 utils::pow_integer<3>
#define pow4 utils::pow_integer<4>

double Waves::compute_constant_CR_source_term_1D() {
  const double x_min = par.source_pmin / cgs::me_c;
  const double x_cut = par.source_cutoff / cgs::me_c;
  const double I = pl_integral(par.source_slope, x_min, x_cut);
  const double R = par.tube_radius;
  double value = initial_luminosity(par.source_luminosity_today, par.source_age, par.source_tdecay);
  value /= 4.0 * pow2(M_PI) * cgs::c_light * pow4(cgs::me_c) * pow2(R) * I;
  return value;
}

double Waves::compute_constant_CR_source_term_3D() {
  const double x_min = par.source_pmin / cgs::me_c;
  const double x_cut = par.source_cutoff / cgs::me_c;
  const double I = pl_integral(par.source_slope, x_min, x_cut);
  double value = initial_luminosity(par.source_luminosity_today, par.source_age, par.source_tdecay);
  value /= 4. * M_PI * cgs::c_light * pow4(cgs::me_c) * I;
  return value;
}

double Waves::compute_constant_CR_source_term() {
  return (par.do_3D) ? compute_constant_CR_source_term_3D() : compute_constant_CR_source_term_1D();
}

void Waves::build_CR_source_term() {
  Q_cr.set_grid_size(par.p_size, par.z_size);
  const double q0 = compute_constant_CR_source_term();
  const double size = par.source_size;
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double F = pl_with_cutoff(p.at(ip), par.source_slope, par.source_pmin, par.source_cutoff);
    for (size_t iz = 0; iz < z.size(); ++iz) {
      const double z_abs = fabs(z.at(iz));
      const double G = (par.do_3D) ? gaussian_3D(z_abs, size) : gaussian_1D(z_abs, size);
      Q_cr.get(ip, iz) = q0 * G * F;
    }
  }
  Q_cr.show_grid("Q_CR", 1. / pow3(cgs::GV) / pow3(cgs::cm) / cgs::second);
}

}  // namespace CRWAVES