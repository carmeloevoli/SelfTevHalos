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

double source_evolution(const double& t_now, const double& t_decay) {
	return std::pow(1. + t_now / t_decay, -2.);
}
