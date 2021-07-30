// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>

double S_i(const double& T_i, const double& gamma_e) {
  double value = 45. * pow2(mks::electron_mass_c2) / 64. / pow2(M_PI) / pow2(mks::k_boltzmann * T_i);
  return value / (value + pow2(gamma_e));
}

void Waves::build_energy_losses() {
  for (size_t ip = 0; ip < p.size(); ++ip) {
    double gamma_e = p.at(ip) / mks::electron_mass_c;
    double energy_density = 0.26 * mks::eV / mks::cm3 * S_i(2.7 * mks::K, gamma_e);  // CMB
    energy_density += 0.30 * mks::eV / mks::cm3 * S_i(20 * mks::K, gamma_e);         // IR
    energy_density += 0.30 * mks::eV / mks::cm3 * S_i(5000 * mks::K, gamma_e);       // star
    // energy_density += 0.10 * eV / cm3 * S_i(20000 * K, gamma_e); // UV
    energy_density += magnetic_energy_density;
    dp_dt.get(ip) = -4. / 3. * mks::sigma_th * energy_density * pow2(gamma_e);
  }
  dp_dt.show_grid("dp_dt", mks::GeV_c / mks::second);
}

}  // namespace CRWAVES