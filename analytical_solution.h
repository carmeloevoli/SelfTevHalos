#ifndef ANALYTICAL_SOLUTION_H_
#define ANALYTICAL_SOLUTION_H_

#include <cmath>
#include "params.h"
#include "units.h"

using namespace std;

class AnalyticSolution {
public:
	AnalyticSolution() {
	}
	~AnalyticSolution() {
	}
	double f(const double& r, const double& t, const double& p) {
		double D_ = D(p);
		double r_ = max(r, 0.01 * pc);
		double A = Q(p) / 4. / M_PI / D_;
		double B = 2. * sqrt(D_ * t);
		return A / r_ * erfc(r / B);
	}
	double df_dz(const double& r, const double& t, const double& p) {
		double D_ = D(p);
		double r_ = max(r, 0.01 * pc);
		double A = Q(p) / 4. / M_PI / D_;
		double B = 2. * sqrt(D_ * t);
		double value = -A / pow2(r_) * erfc(r / B) - 2. * A / r_ / B / sqrt(M_PI) * exp(-pow2(r / B));
		return fabs(value);
	}
	void set_Q0(const double& Q0) {
		this->Q0 = Q0;
	}
	void set_D0(const double& D0) {
		this->D0 = D0;
	}
	void set_D0_ref(const double& D0_ref) {
		this->D0_ref = D0_ref;
	}
	void set_delta(const double& delta) {
		this->delta = delta;
	}
	void set_alpha(const double& alpha) {
		this->alpha = alpha;
	}
	void set_p_cutoff(const double& p_cutoff) {
		this->p_cutoff = p_cutoff;
	}
	void set_params(const Params& params, const double& Q0) {
		set_D0(params.D_gal());
		set_D0_ref(params.D_gal_ref());
		set_delta(1. / 3.);
		set_alpha(params.alpha());
		set_p_cutoff(params.source_cutoff());
		set_Q0(Q0);
	}
private:
	double Q0 = 0;
	double D0 = 3e28 * cm2 / s;
	double D0_ref = GeV_c;
	double delta = 1. / 3.;
	double alpha = 3.5;
	double p_cutoff = 100 * TeV_c;

	double Q(const double& p) {
		return Q0 * pow(p / electron_mass_c, -alpha) * exp(-pow2(p / p_cutoff));
	}
	double D(const double& p) {
		return D0 * pow(p / D0_ref, delta);
	}
};

#endif /* ANALYTICAL_SOLUTION_H_ */
