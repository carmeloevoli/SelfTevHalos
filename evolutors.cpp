#include "waves.h"
#include "omp.h"
#define OMP_NUM_THREADS 4

double Waves::source_evolution(const double& t_now) {
	return pow(1. + t_now / params.source_tdecay(), -2.);
}

void Waves::evolve_f_diffusive(const size_t& number_of_operators, const double& t_now) {

	double dt_half = 0.5 * dt;
	double dz = fabs(z.at(1) - z.at(0));
	double vA_abs_dz = abs(params.vA_infty()) / abs(z.at(1) - z.at(0));

#pragma omp parallel for
	for (int ip = 0; ip < p_size - 1; ++ip) {
		vector<double> rhs(z_size - 2);
		vector<double> central_diagonal(z_size - 2);
		vector<double> upper_diagonal(z_size - 3);
		vector<double> lower_diagonal(z_size - 3);
		vector<double> fcr_up(z_size - 2);

		for (int iz = 1; iz < z_size - 1; ++iz) {
			double Dz = D_zz.get(ip, iz);
			double DzUp = (iz < z_size - 1) ? D_zz.get(ip, iz + 1) : Dz;
			double DzDo = (iz > 0) ? D_zz.get(ip, iz - 1) : Dz;

			double Dz_dz2 = Dz / pow2(dz);
			double dDz_4_dz2 = (DzUp - DzDo) / 4. / pow2(dz);
			double r = max(fabs(z.at(iz)), .2 * pc);
			double Dz_dz_r = Dz / dz / r;

			double UZ = Dz_dz2 + dDz_4_dz2; // + Dz_dz_r;
			double CZ = 2. * Dz_dz2;
			double LZ = Dz_dz2 - dDz_4_dz2; // - Dz_dz_r;

			/* advective term */
			//UZ += 0;
			//CZ += vA_abs_dz;
			//LZ += vA_abs_dz;
			//} else if (iz < z.get_idx()) {
			//	LZ += 0;
			//	CZ += -vA_abs_dz;
			//	UZ += vA_abs_dz;
			//} else {
			//	LZ += -vA_abs_dz / 2.;
			//	CZ += -vA_abs_dz;
			//	UZ += -vA_abs_dz / 2.;
			//}*/
			/* end advective term */

			central_diagonal.at(iz - 1) = 1. + dt_half * CZ;
			if (iz != z_size - 2) {
				upper_diagonal.at(iz - 1) = -dt_half * UZ;
			}
			/*if (iz == 0) {
				upper_diagonal.at(iz) += -dt_half * LZ;
			}*/
			if (iz != 1) {
				lower_diagonal.at(iz - 2) = -dt_half * LZ;
			}
			rhs.at(iz - 1) = fcr.get(ip, iz) * (2. - central_diagonal.at(iz - 1));
			rhs.at(iz - 1) += dt * Q_cr.get(ip, iz) * source_evolution(t_now) / (double) number_of_operators;

			if (iz != z_size - 2) {
				rhs.at(iz - 1) -= fcr.get(ip, iz + 1) * upper_diagonal.at(iz - 1);
			}
			if (iz != 1) {
				rhs.at(iz - 1) -= fcr.get(ip, iz - 1) * lower_diagonal.at(iz - 2);
			}
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (int iz = 1; iz < z_size - 1; ++iz) {
			double value = .5 * (fcr_up.at(iz - 1) + fcr_up.at(z_size - 2 - iz));
			fcr.get(ip, iz) = max(value, 0.);
		}
	} // for
}

void Waves::evolve_f_adiabatic(const size_t& number_of_operators, const double& t_now) {

	double dt_half = 0.5 * dt;
	double dvadz = 0; // (v_A.get(iz + 1) - v_A.get(iz - 1)) / (z.at(iz + 1) - z.at(iz - 1));

#pragma omp parallel for
	for (int iz = 1; iz < z_size - 1; ++iz) {
		vector<double> rhs(p_size - 1);
		vector<double> central_diagonal(p_size - 1);
		vector<double> upper_diagonal(p_size - 2);
		vector<double> lower_diagonal(p_size - 2);
		vector<double> fcr_up(p_size - 1);

		for (int ip = 0; ip < p_size - 1; ++ip) {

			double b_i = dvadz / 3.0 / (p.at(ip + 1) / p.at(ip) - 1.0);
			double L = 0, C = 0, U = 0;

			C += b_i;
			U += b_i;

			C -= dpdt.get(ip, iz) / (p.at(ip + 1) - p.at(ip));
			U -= pow2(p.at(ip + 1) / p.at(ip)) * dpdt.get(ip + 1, iz) / (p.at(ip + 1) - p.at(ip));

			central_diagonal.at(ip) = 1. + dt_half * C;
			if (ip != 0) {
				lower_diagonal.at(ip - 1) = -dt_half * L;
			}
			if (ip != p_size - 2) {
				upper_diagonal.at(ip) = -dt_half * U;
			}
			rhs.at(ip) = fcr.get(ip, iz) * (2. - central_diagonal.at(ip));
			rhs.at(ip) += dt * Q_cr.get(ip, iz) * source_evolution(t_now) / (double) number_of_operators;
			if (ip != p_size - 2) {
				rhs.at(ip) -= fcr.get(ip + 1, iz) * upper_diagonal.at(ip);
			}
			if (ip != 0) {
				rhs.at(ip) -= fcr.get(ip - 1, iz) * lower_diagonal.at(ip - 1);
			}

		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (int ip = 0; ip < p_size - 1; ++ip) {
			double value = fcr_up.at(ip);
			fcr.get(ip, iz) = max(value, 0.);
		}
	} // for
}

void Waves::evolve_f_adiabatic_2(const size_t& number_of_operators, const double& t_now) {

	double dt_half = 0.5 * dt;
	double dvadz = 0; // (v_A.get(iz + 1) - v_A.get(iz - 1)) / (z.at(iz + 1) - z.at(iz - 1));

	for (int iz = 0; iz < z_size - 1; ++iz) {
		vector<double> rhs(p_size - 1);
		vector<double> central_diagonal(p_size - 1);
		vector<double> upper_diagonal(p_size - 2);
		vector<double> lower_diagonal(p_size - 2);
		vector<double> fcr_up(p_size - 1);

		for (int ip = 1; ip < p_size; ++ip) {

			double b_i = 0; // dvadz / 3.0 / (p.at(ip + 1) / p.at(ip) - 1.0);
			double L = 0, C = 0, U = 0;

			C += b_i;
			U += b_i;

			C -= dpdt.get(ip, iz) / (p.at(ip) - p.at(ip - 1));
			L -= pow2(p.at(ip - 1) / p.at(ip)) * dpdt.get(ip - 1, iz) / (p.at(ip) - p.at(ip - 1));

			central_diagonal.at(ip - 1) = 1. + dt_half * C;
			if (ip != 1) {
				lower_diagonal.at(ip - 2) = -dt_half * L;
			}
			if (ip != p_size - 1) {
				upper_diagonal.at(ip - 1) = -dt_half * U;
			}
			rhs.at(ip - 1) = fcr.get(ip, iz) * (2. - central_diagonal.at(ip - 1));
			rhs.at(ip - 1) += dt * Q_cr.get(ip, iz) * source_evolution(t_now) / (double) number_of_operators;
			if (ip != p_size - 1) {
				rhs.at(ip - 1) -= fcr.get(ip + 1, iz) * upper_diagonal.at(ip - 1);
			}
			if (ip != 1) {
				rhs.at(ip - 1) -= fcr.get(ip - 1, iz) * lower_diagonal.at(ip - 2);
			}
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (int ip = 1; ip < p_size; ++ip) {
			double value = fcr_up.at(ip - 1);
			fcr.get(ip, iz) = max(value, 0.);
		}
	} // for
}

void Waves::evolve_f_adiabatic_2nd(const size_t& number_of_operators, const double& t_now) {
	for (int iz = 0; iz < z_size - 1; ++iz) {
		vector<double> rhs(p_size - 1);
		vector<double> central_diagonal(p_size - 1);
		vector<double> upper_diagonal(p_size - 2);
		vector<double> lower_diagonal(p_size - 2);
		vector<double> fcr_up(p_size - 1);

		for (int ip = 0; ip < p_size - 1; ++ip) {

			double C = -dpdt.get(ip, iz) / (p.at(ip + 1) - p.at(ip));
			double U = -dpdt.get(ip + 1, iz) / (p.at(ip + 1) - p.at(ip));

			central_diagonal.at(ip) = 1. + dt * C;
			if (ip != p_size - 2) {
				upper_diagonal.at(ip) = 1. - dt * U;
			}
			if (ip != 0) {
				lower_diagonal.at(ip - 1) = 0.;
			}
			rhs.at(ip) = (1. - dt * C) * fcr.get(ip, iz) + (1. + dt * U) * fcr.get(ip + 1, iz)
					+ dt * 0.5 * (Q_cr.get(ip + 1, iz) + Q_cr.get(ip, iz)) / (double) number_of_operators;
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (int ip = 0; ip < p_size - 1; ++ip) {
			double value = fcr_up.at(ip);
			fcr.get(ip, iz) = max(value, 0.);
		}
	} // for iz
}

void Waves::evolve_waves_advectice(const size_t& number_of_operators) {
	double vA_abs_dz = abs(params.vA_infty()) / abs(z.at(1) - z.at(0));
	double dt_half = 0.5 * dt;
	double factor_damping = pow(2. * params.ck(), -1.5) * params.vA_infty();
	double factor_growth = 2. * M_PI / 3. * c_light * params.vA_infty() / params.magnetic_energy_density();

#pragma omp parallel for
	for (size_t ip = 0; ip < p_size - 1; ++ip) {
		vector<double> rhs(z_size - 1);
		vector<double> central_diagonal(z_size - 1);
		vector<double> upper_diagonal(z_size - 2);
		vector<double> lower_diagonal(z_size - 2);
		vector<double> W_up(z_size - 1);

		double k_ = 1. / larmor_radius(p.at(ip), params.magnetic_field());

		for (size_t iz = 0; iz < z_size - 1; ++iz) {
			double L, C, U;
			L = vA_abs_dz;
			C = vA_abs_dz;
			U = 0.;

			central_diagonal.at(iz) = 1. + dt_half * C;
			if (iz != 0) {
				lower_diagonal.at(iz - 1) = -dt_half * L;
			}
			if (iz != z_size - 2) {
				upper_diagonal.at(iz) = -dt_half * U;
			}
			if (iz == 0) {
				upper_diagonal.at(iz) += -dt_half * L;
			}

			rhs.at(iz) = W_sg.get(ip, iz) * (2. - central_diagonal.at(iz));
			rhs.at(iz) += (iz != 0) ? dt_half * L * W_sg.get(ip, iz - 1) : 0.;
			rhs.at(iz) += (iz != z_size - 2) ? dt_half * U * W_sg.get(ip, iz + 1) : 0;

			double Gamma_D = factor_damping * pow(k_, 1.5) * sqrt(W_sg.get(ip, iz));
			double Gamma_D_gal = factor_damping * pow(k_, 1.5) * sqrt(W_ext.get(ip, iz));
			double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * dfdz.get(ip, iz);

			double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ext.get(ip, iz);

			rhs.at(iz) += dt * Q_w / (double) number_of_operators;
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, W_up);

		for (size_t iz = 0; iz < z_size - 1; ++iz) {
			double value = W_up.at(iz); // min(W_up.at(iz), 1. / k_);
			W_sg.get(ip, iz) = max(value, 0.);
		}
	}
}

void Waves::evolve_waves_noadvectice(const size_t& number_of_operators) {
	double factor_damping = pow(2. * params.ck(), -1.5) * params.vA_infty();
	double factor_growth = 2. * M_PI / 3. * c_light * params.vA_infty() / params.magnetic_energy_density();

#pragma omp parallel for
	for (size_t ip = 0; ip < p_size; ++ip) {
		double k_ = 1. / larmor_radius(p.at(ip), params.magnetic_field());
		for (size_t iz = 0; iz < z_size; ++iz) {
			double Gamma_D = factor_damping * pow(k_, 1.5) * sqrt(W_sg.get(ip, iz));
			double Gamma_D_gal = factor_damping * pow(k_, 1.5) * sqrt(W_ext.get(ip, iz));
			double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * dfdz.get(ip, iz);
			double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ext.get(ip, iz);
			double value = W_sg.get(ip, iz) + dt * Q_w;
			W_sg.get(ip, iz) = max(value, 0.);
		}
	}
}
