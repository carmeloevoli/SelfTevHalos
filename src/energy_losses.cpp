#include "waves.h"

double S_i(const double& T_i, const double& gamma_e) {
	double value = 45. * pow2(electron_mass_c2) / 64. / pow2(M_PI) / pow2(k_boltzmann * T_i);
	return value / (value + pow2(gamma_e));
}

void Waves::build_energy_losses() {
	dp_dt.set_grid_size(p.get_size(), 1);
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		double gamma_e = p.at(ip) / electron_mass_c;
		double energy_density = 0.26 * eV / cm3 * S_i(2.7 * K, gamma_e); // CMB
		energy_density += 0.30 * eV / cm3 * S_i(20 * K, gamma_e); // IR
		energy_density += 0.30 * eV / cm3 * S_i(5000 * K, gamma_e); // star
		//energy_density += 0.10 * eV / cm3 * S_i(20000 * K, gamma_e); // UV
		energy_density += magnetic_energy_density;
		dp_dt.get(ip) = -4. / 3. * sigma_th * energy_density * pow2(gamma_e);
	}
	dp_dt.show_grid("dp_dt", 1.);
}
