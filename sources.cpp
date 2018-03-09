#include "waves.h"

double f(double x, void * params) {
	double alpha = *(double *) params;
	double x_cutoff = *((double *) params + 1);
	double f = pow(x, 3. - alpha) * exp(-pow2(x / x_cutoff));
	return f;
}

double I_of_alpha(double alpha, double x_min, double x_cutoff) {
	int LIMIT = 10000;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
	double result, error;
	double params[2] = { alpha, x_cutoff };
	gsl_function F;
	F.function = &f;
	F.params = &params;
	gsl_integration_qag(&F, x_min, 2.0 * x_cutoff, 0, 1e-3, LIMIT, 3, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}

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
	double R = 1 * pc; // TODO move to params
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
