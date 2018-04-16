#include "waves.h"

#define sgn(A) ((A < 0.) ? -1. : 1.)

Waves::Waves(const Params& par_) {
	par = par_;
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
	if (par.do_3D.get())
		z.build_lin_axis(0, halo_size, z_size);
	else
		z.build_lin_axis(-halo_size, halo_size, z_size);
	this->z_size = z.get_size();
	z.set_reference_value(0);
	z.show_axis("z", kpc);
	dz = fabs(z.at(1) - z.at(0));
}

void Waves::build_W_ISM() {
	W_ISM.set_grid_size(p.get_size(), z.get_size());
	double k_norm = 1. / larmor_radius(par.D_gal_ref.get(), par.magnetic_field.get());
	double eta_B = c_light / 2. / k_norm / par.D_gal.get() * pow(k_norm / par.k0.get(), 2. / 3.);
	std::cout << " - eta_B = " << eta_B << "\n";
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		for (size_t iz = 0; iz < z.get_size(); ++iz) {
			double value = 2.0 * eta_B / 3.0 / par.k0.get();
			double k = 1. / larmor_radius(p.at(ip), par.magnetic_field.get());
			if (k > par.k0.get())
				value *= pow(k / par.k0.get(), -5. / 3.);
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
		//v_A.get(j) = sgn(z.at(j)) * par.vA_infty();
		v_A.get(j) = std::tanh(z.at(j) / pc) * par.vA_infty.get();
	}
	v_A.show_grid("v_A", km / s);
}

void Waves::build_analytical_solution() {
	solution.set_params(par, compute_constant_CR_source_term());
}
