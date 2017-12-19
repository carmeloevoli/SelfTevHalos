#include "waves.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

void Waves::compute_D_zz_test() {
	double H = z.get_max();
	double D0 = pc * pc / year;
	for (size_t ip = 0; ip < p_size; ++ip) {
		for (size_t iz = 0; iz < z_size; ++iz) {
			double r = max(fabs(z.at(iz)), 5. * pc);
			D_zz.get(ip, iz) = D0 * (1. + pow2(r / H));
		}
	}
}

void Waves::compute_D_zz() {
	for (size_t ip = 0; ip < p_size; ++ip) {
		vector<double> D_temp(z_size);
		for (size_t iz = 1; iz < z_size - 1; ++iz) {
			double k = 1. / larmor_radius(p.at(ip), params.magnetic_field());
			double value = beta(p.at(ip)) * c_light / 3.0 / pow2(k) / W_sg.get(ip, iz);
			D_temp.at(iz) = value;
		}
		for (size_t iz = 1; iz < z_size - 1; ++iz) {
			D_zz.get(ip, iz) = (D_temp.at(iz - 1) + D_temp.at(iz) + D_temp.at(iz + 1)) / 3.0;
		}
		D_zz.get(ip, 0) = D_zz.get(ip, 1);
		D_zz.get(ip, z_size - 1) = D_zz.get(ip, z_size - 2);
	}
}

void Waves::smooth_D_zz(const int& smoothing_radius) {
	//std::cout << "smoothing Dzz..." << "\n";
	double dz = z.at(1) - z.at(0);
	double sigma2 = pow2(smoothing_radius * dz);
	TGrid2D<double> D_zz_smoothed;
	D_zz_smoothed.set_grid_size(p_size, z_size);
	for (size_t ip = 0; ip < p_size; ++ip) {
		for (size_t iz = 0; iz < z_size; ++iz) {
			double D_smoothed = 0;
			for (size_t iz_prime = 0; iz_prime < z_size; ++iz_prime) {
				double G_filter = 1. / sqrt(2. * M_PI * sigma2) * exp(-pow2(z.at(iz_prime) - z.at(iz)) / 2. / sigma2);
				D_smoothed += dz * G_filter * D_zz.get(ip, iz_prime);
			}
			D_zz_smoothed.get(ip, iz) = D_smoothed;
		}
	}
	for (size_t ip = 0; ip < p_size; ++ip) {
		for (size_t iz = 0; iz < z_size; ++iz) {
			D_zz.get(ip, iz) = D_zz_smoothed.get(ip, iz);
		}
	}
}

void Waves::smooth_W(const int& smoothing_radius) {
	std::cout << "smoothing Waves..." << "\n";
	double dz = z.at(1) - z.at(0);
	double sigma2 = pow2(smoothing_radius * dz);
	TGrid2D<double> W_smoothed;
	W_smoothed.set_grid_size(p_size, z_size);
	for (size_t ip = 0; ip < p_size; ++ip) {
		for (size_t iz = 0; iz < z_size; ++iz) {
			double W_smoothed_ = 0;
			for (size_t iz_prime = iz - 10; iz_prime <= iz + 10; ++iz_prime) {
				if (iz_prime > -1 && iz_prime < z_size) {
					double G_filter = 1. / sqrt(2. * M_PI * sigma2) * exp(-pow2(z.at(iz_prime) - z.at(iz)) / 2. / sigma2);
					W_smoothed_ += dz * G_filter * W_sg.get(ip, iz_prime);
				}
			}
			W_smoothed.get(ip, iz) = W_smoothed_;
		}
	}
	for (size_t ip = 0; ip < p_size; ++ip) {
		for (size_t iz = 0; iz < z_size; ++iz) {
			W_sg.get(ip, iz) = W_smoothed.get(ip, iz);
		}
	}
}

double Waves::compute_dfdz() {
	size_t nz = z.get_size();
	size_t np = p.get_size();
	double value = 0;
	double dz = fabs(z.at(1) - z.at(0));
	for (size_t iz = 1; iz < nz - 2; ++iz) {
		for (size_t ip = 0; ip < np; ++ip) {
			//dfdz.get(ip, iz) = fcr.get(ip, iz) / fabs(z.at(iz) + pc); // max(abs(value), 1e-100);
			//dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
			if (iz == z.get_idx()) {
				value = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
			} else {
				value = fabs((fcr.get(ip, iz + 1) - fcr.get(ip, iz - 1)) / 2. / dz);
				//dfdz.get(ip, iz) = fabs(-fcr.get(ip, iz + 2) + 4.0 * fcr.get(ip, iz + 1) - 3.0 * fcr.get(ip, iz)) / 2. / dz;
				//dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
			}
			dfdz.get(ip, iz) = max(value, 0.);
		}

	}
	return 0;
}
