#include "waves.h"

#define sgn(A) ((A < 0.) ? -1. : 1.)

Waves::Waves() {
	par.print();
}

Waves::~Waves() {
}

void Waves::build_p_axis(const double& p_min, const double& p_max, const size_t& p_size) {
	p.build_log_axis(p_min, p_max, p_size);
	this->p_size = p.get_size();
	p.set_reference_value(sqrt(p_min * p_max));
	p.show_axis("p", GeV_c);
	dlnp = std::log(p.at(1) / p.at(0));
}

void Waves::build_z_axis(const double& halo_size, const size_t& z_size) {
	z.build_lin_axis(-halo_size, halo_size, z_size);
	this->z_size = z.get_size();
	z.set_reference_value(0);
	z.show_axis("z", kpc);
	dz = fabs(z.at(1) - z.at(0));
}

double Waves::compute_constant_CR_source_term_1D() {
	double I = I_of_alpha(par.alpha(), par.source_pmin() / electron_mass_c, par.source_cutoff() / electron_mass_c);
	double out = par.spin_down_luminosity() * pow2(1. + par.age() / par.source_tdecay());
	double R = 1. * pc; // TODO move to params
	out /= 4.0 * pow2(M_PI) * c_light * pow4(electron_mass_c) * pow2(R) * I;
	return out;
}

double Waves::compute_constant_CR_source_term_3D() {
	double I = I_of_alpha(par.alpha(), par.source_pmin() / electron_mass_c, par.source_cutoff() / electron_mass_c);
	double out = par.spin_down_luminosity() * pow2(1. + par.age() / par.source_tdecay());
	out /= 4. * M_PI * c_light * pow4(electron_mass_c) * I;
	return out;
}

void Waves::build_CR_source_term() {
	Q_cr.set_grid_size(p.get_size(), z.get_size());
	double q0 = (par.do_3D()) ? compute_constant_CR_source_term_3D() : compute_constant_CR_source_term_1D();
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

void Waves::build_W_ISM() {
	W_ISM.set_grid_size(p.get_size(), z.get_size());
	double k_norm = 1. / larmor_radius(par.D_gal_ref(), par.magnetic_field());
	double eta_B = c_light / 2. / k_norm / par.D_gal() * pow(k_norm / par.k0(), 2. / 3.);
	std::cout << " - eta_B = " << eta_B << "\n";
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		for (size_t iz = 0; iz < z.get_size(); ++iz) {
			double value = 2.0 * eta_B / 3.0 / par.k0();
			double k = 1. / larmor_radius(p.at(ip), par.magnetic_field());
			if (k > par.k0())
				value *= pow(k / par.k0(), -5. / 3.);
			W_ISM.get(ip, iz) = value;
		}
	}
	W_ISM.show_grid("W_ISM", 1.);
}

void Waves::build_W_sg() {
	W_sg.set_grid_size(p.get_size(), z.get_size());
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		for (size_t iz = 0; iz < z.get_size(); ++iz) {
			W_sg.get(ip, iz) = W_ISM.get(ip, iz);
		}
	}
	W_sg.show_grid("W_sg", 1.);
}

void Waves::build_f_cr() {
	df_dz.set_grid_size(p.get_size(), z.get_size());
	f_cr.set_grid_size(p.get_size(), z.get_size());
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		for (size_t iz = 0; iz < z.get_size(); ++iz) {
			df_dz.get(ip, iz) = 0;
			f_cr.get(ip, iz) = 0;
		}
	}
	df_dz.show_grid("dfdz", 1.);
	f_cr.show_grid("fp", 1.);
}

void Waves::build_D_zz() {
	D_zz.set_grid_size(p.get_size(), z.get_size());
	compute_D_zz();
	D_zz.show_grid("D_zz", 1.);
}

void Waves::build_v_A() {
	v_A.set_grid_size(1, z.get_size());
	for (size_t j = 0; j < z.get_size(); ++j) {
		v_A.get(j) = sgn(z.at(j)) * par.vA_infty();
	}
	v_A.show_grid("v_A", km / s);
}

void Waves::build_energy_losses() {
	dp_dt.set_grid_size(p.get_size(), 1);
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		double gamma_e = p.at(ip) / electron_mass_c;
		double energy_density = 0.26 * eV / cm3 * S_i(2.7 * K, gamma_e); // CMB
		energy_density += 0.30 * eV / cm3 * S_i(20 * K, gamma_e); // IR
		energy_density += 0.30 * eV / cm3 * S_i(5000 * K, gamma_e); // star
		//energy_density += 0.10 * eV / cm3 * S_i(20000 * K, gamma_e); // UV
		energy_density += par.magnetic_energy_density();
		dp_dt.get(ip) = -4. / 3. * sigma_th * energy_density * pow2(gamma_e);
	}
	dp_dt.show_grid("dp_dt", 1.);
}
