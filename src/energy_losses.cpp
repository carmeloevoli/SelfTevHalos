// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>

double S_i(const double& T_i, const double& gamma_e) {
  double value = 45. * pow2(cgs::me_c2) / 64. / pow2(M_PI) / pow2(cgs::k_B * T_i);
  return value / (value + pow2(gamma_e));
}

void Waves::build_energy_losses() {
  dp_dt.resize(par.p_size);
  for (size_t ip = 0; ip < p.size(); ++ip) {
    double gamma_e = p.at(ip) / cgs::me_c;
    double energy_density = 0.26 * cgs::eV / cgs::cm3 * S_i(2.7 * cgs::K, gamma_e);  // CMB
    energy_density += 0.30 * cgs::eV / cgs::cm3 * S_i(20 * cgs::K, gamma_e);         // IR
    energy_density += 0.30 * cgs::eV / cgs::cm3 * S_i(5000 * cgs::K, gamma_e);       // star
    // energy_density += 0.10 * eV / cm3 * S_i(20000 * K, gamma_e); // UV
    energy_density += magnetic_energy_density(par.magnetic_field);
    dp_dt.at(ip) = -4. / 3. * cgs::sigma_th * energy_density * pow2(gamma_e);
  }
  utils::infoAxis(dp_dt, "dp_dt", cgs::GV / cgs::second);
}

}  // namespace CRWAVES