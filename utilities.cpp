#include "utilities.h"

size_t find_nearest(const std::vector<double>& x, const double& x0) {
	double min_distance = 1e100;
	size_t min_value = 0;
	for (size_t i = 0; i < x.size(); ++i)
		if (fabs(x[i] - x0) < min_distance) {
			min_distance = fabs(x[i] - x0);
			min_value = i;
		}
	return min_value;
}

void print_vector(const std::vector<double>& x) {
	for (size_t i = 0; i < x.size(); ++i)
		std::cout << i << "\t" << x.at(i) << "\n";
}

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

double beta(const double& p) {
	double gamma2 = 1. + pow2(p / electron_mass_c);
	return sqrt(1. - 1. / gamma2);
}

double larmor_radius(const double& p, const double& B) {
	return p / elementary_charge / B;
}

double inverse_larmor_radius(const double& k, const double& B) {
	return elementary_charge * B / k;
}

double spectrum(const double& p, const double& alpha, const double& p_min, const double& p_cutoff) {
	return (p > p_min) ? pow(p / electron_mass_c, -alpha) * exp(-pow2(p / p_cutoff)) : 0.;
}

double source_profile_1D(const double& z, const double& size) {
	return pow(2.0 * M_PI * pow2(size), -1. / 2.) * exp(-0.5 * pow2(z / size));
}

double source_profile_3D(const double& z, const double& size) {
	return pow(2.0 * M_PI * pow2(size), -3. / 2.) * exp(-0.5 * pow2(z / size));
}

double S_i(const double& T_i, const double& gamma_e) {
	double value = 45. * pow2(electron_mass_c2) / 64. / pow2(M_PI) / pow2(k_boltzmann * T_i);
	return value / (value + pow2(gamma_e));
}

double source_evolution(const double& t_now, const double& t_decay) {
	return std::pow(1. + t_now / t_decay, -2.);
}
