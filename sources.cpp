#include "waves.h"

double source_profile_1D(const double& z, const double& size) {
	return pow(2.0 * M_PI * pow2(size), -1. / 2.) * exp(-0.5 * pow2(z / size));
}

double source_profile_3D(const double& z, const double& size) {
	return pow(2.0 * M_PI * pow2(size), -3. / 2.) * exp(-0.5 * pow2(z / size));
}

double spectrum(const double& p, const double& alpha, const double& p_min, const double& p_cutoff) {
	return (p > p_min) ? pow(p / electron_mass_c, -alpha) * exp(-pow2(p / p_cutoff)) : 0.;
}

double Waves::compute_constant_CR_source_term_1D() {
	double I = I_of_alpha(par.alpha(), par.source_pmin() / electron_mass_c, par.source_cutoff() / electron_mass_c);
	double out = par.spin_down_luminosity() * pow2(1. + par.age() / par.source_tdecay());
	double R = 0.5 * pc; // TODO move to params
	out /= 4.0 * pow2(M_PI) * c_light * pow4(electron_mass_c) * pow2(R) * I;
	return out;
}

double Waves::compute_constant_CR_source_term_3D() {
	double I = I_of_alpha(par.alpha(), par.source_pmin() / electron_mass_c, par.source_cutoff() / electron_mass_c);
	double out = par.spin_down_luminosity() * pow2(1. + par.age() / par.source_tdecay());
	out /= 4. * M_PI * c_light * pow4(electron_mass_c) * I;
	return out;
}

double Waves::compute_constant_CR_source_term() {
	return (par.do_3D()) ? compute_constant_CR_source_term_3D() : compute_constant_CR_source_term_1D();
}

void Waves::build_CR_source_term() {
	Q_cr.set_grid_size(p.get_size(), z.get_size());
	double q0 = compute_constant_CR_source_term();
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		double F = spectrum(p.at(ip), par.alpha(), par.source_pmin(), par.source_cutoff());
		for (size_t iz = 0; iz < z.get_size(); ++iz) {
			double z_ = fabs(z.at(iz));
			double G = (par.do_3D()) ? source_profile_3D(z_, par.source_size()) : source_profile_1D(z_, par.source_size());
			Q_cr.get(ip, iz) = q0 * G * F;
		}
	}
	Q_cr.show_grid("Q_CR", 1.);
}



