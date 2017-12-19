#include "waves.h"
#include "omp.h"

#define OMP_NUM_THREADS 4

void Waves::evolve_f_in_z(const size_t& number_of_operators, const double& t_now) {
	double vA_abs_dz = abs(par.vA_infty()) / abs(z.at(1) - z.at(0));

#pragma omp parallel for
	for (int ip = 0; ip < p_size - 1; ++ip) {
		std::vector<double> rhs(z_size - 2);
		std::vector<double> central_diagonal(z_size - 2);
		std::vector<double> upper_diagonal(z_size - 3);
		std::vector<double> lower_diagonal(z_size - 3);
		std::vector<double> fcr_up(z_size - 2);

		for (int iz = 1; iz < z_size - 1; ++iz) {
			double Dz = D_zz.get(ip, iz);
			double DzUp = (iz < z_size - 1) ? D_zz.get(ip, iz + 1) : Dz;
			double DzDo = (iz > 0) ? D_zz.get(ip, iz - 1) : Dz;

			double Dz_dz2 = Dz / pow2(dz);
			double dDz_4_dz2 = (DzUp - DzDo) / 4. / pow2(dz);
			double r = std::max(fabs(z.at(iz)), .2 * pc); // TODO smaller radius
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
			rhs.at(iz - 1) = f_cr.get(ip, iz) * (2. - central_diagonal.at(iz - 1));
			if (iz != z_size - 2) {
				rhs.at(iz - 1) -= f_cr.get(ip, iz + 1) * upper_diagonal.at(iz - 1);
			}
			if (iz != 1) {
				rhs.at(iz - 1) -= f_cr.get(ip, iz - 1) * lower_diagonal.at(iz - 2);
			}
			double Q = Q_cr.get(ip, iz) * source_evolution(t_now, par.source_tdecay());
			rhs.at(iz - 1) += dt * Q / (double) number_of_operators;
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (int iz = 1; iz < z_size - 1; ++iz) {
			double value = .5 * (fcr_up.at(iz - 1) + fcr_up.at(z_size - 2 - iz));
			f_cr.get(ip, iz) = std::max(value, 0.);
		}
	} // for
}

void Waves::evolve_f_in_p(const size_t& number_of_operators, const double& t_now) {

#pragma omp parallel for
	for (int iz = 1; iz < z_size - 1; ++iz) {
		std::vector<double> rhs(p_size - 1);
		std::vector<double> central_diagonal(p_size - 1);
		std::vector<double> upper_diagonal(p_size - 2);
		std::vector<double> lower_diagonal(p_size - 2);
		std::vector<double> fcr_up(p_size - 1);
		double dvadz = (v_A.get(iz + 1) - v_A.get(iz - 1)) / (z.at(iz + 1) - z.at(iz - 1));

		for (int ip = 0; ip < p_size - 1; ++ip) {
			double b_i = 0; // dvadz / 3.0 / (p.at(ip + 1) / p.at(ip) - 1.0);
			double L = 0, C = 0, U = 0;

			C += b_i;
			U += b_i;

			C -= dp_dt.get(ip) / (p.at(ip + 1) - p.at(ip));
			U -= pow2(p.at(ip + 1) / p.at(ip)) * dp_dt.get(ip + 1) / (p.at(ip + 1) - p.at(ip));

			central_diagonal.at(ip) = 1. + dt_half * C;
			if (ip != 0) {
				lower_diagonal.at(ip - 1) = -dt_half * L;
			}
			if (ip != p_size - 2) {
				upper_diagonal.at(ip) = -dt_half * U;
			}
			rhs.at(ip) = f_cr.get(ip, iz) * (2. - central_diagonal.at(ip));
			if (ip != p_size - 2) {
				rhs.at(ip) -= f_cr.get(ip + 1, iz) * upper_diagonal.at(ip);
			}
			if (ip != 0) {
				rhs.at(ip) -= f_cr.get(ip - 1, iz) * lower_diagonal.at(ip - 1);
			}
			double Q = Q_cr.get(ip, iz) * source_evolution(t_now, par.source_tdecay());
			rhs.at(ip) += dt * Q / (double) number_of_operators;
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (int ip = 0; ip < p_size - 1; ++ip) {
			double value = fcr_up.at(ip);
			f_cr.get(ip, iz) = std::max(value, 0.);
		}
	} // for
}

void Waves::evolve_waves_in_z(const size_t& number_of_operators) {
	double vA_abs_dz = abs(par.vA_infty()) / abs(z.at(1) - z.at(0));

#pragma omp parallel for
	for (size_t ip = 0; ip < p_size - 1; ++ip) {
		std::vector<double> rhs(z_size - 1);
		std::vector<double> central_diagonal(z_size - 1);
		std::vector<double> upper_diagonal(z_size - 2);
		std::vector<double> lower_diagonal(z_size - 2);
		std::vector<double> W_up(z_size - 1);

		double k_ = 1. / larmor_radius(p.at(ip), par.magnetic_field());

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
			double Gamma_D_gal = factor_damping * pow(k_, 1.5) * sqrt(W_ISM.get(ip, iz));
			double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * df_dz.get(ip, iz);
			double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ISM.get(ip, iz);
			rhs.at(iz) += dt * Q_w / (double) number_of_operators;
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, W_up);

		for (size_t iz = 0; iz < z_size - 1; ++iz) {
			double value = W_up.at(iz); // min(W_up.at(iz), 1. / k_);
			W_sg.get(ip, iz) = std::max(value, 0.);
		}
	}
}

void Waves::evolve_waves(const size_t& number_of_operators) {
#pragma omp parallel for
	for (size_t ip = 0; ip < p_size; ++ip) {
		double k_ = 1. / larmor_radius(p.at(ip), par.magnetic_field());
		for (size_t iz = 0; iz < z_size; ++iz) {
			double Gamma_D = factor_damping * pow(k_, 1.5) * sqrt(W_sg.get(ip, iz));
			double Gamma_D_gal = factor_damping * pow(k_, 1.5) * sqrt(W_ISM.get(ip, iz));
			double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * df_dz.get(ip, iz);
			double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ISM.get(ip, iz);
			double value = W_sg.get(ip, iz) + dt * Q_w;
			W_sg.get(ip, iz) = std::max(value, 0.);
		}
	}
}
