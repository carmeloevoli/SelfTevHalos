#include "waves.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

void Waves::compute_D_zz() {
	for (size_t ip = 0; ip < p_size; ++ip) {
		double k = 1. / larmor_radius(p.at(ip), par.magnetic_field());
		double factor = beta(p.at(ip)) * c_light / 3.0 / pow2(k);
		for (size_t iz = 1; iz < z_size - 1; ++iz) {
			double value = factor / W_sg.get(ip, iz);
			D_zz.get(ip, iz) = value;
		}
		D_zz.get(ip, 0) = D_zz.get(ip, 1);
		D_zz.get(ip, z_size - 1) = D_zz.get(ip, z_size - 2);
	}
}

void Waves::compute_dfdz() {
	for (size_t ip = 0; ip < p_size; ++ip) {
		double value = 0;
		for (size_t iz = 1; iz < z_size - 1; ++iz) {
			if (iz == z.get_idx())
				value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
			else
				value = fabs((f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1)) / 2. / dz);
			df_dz.get(ip, iz) = fabs(value);
		}
		df_dz.get(ip, 0) = df_dz.get(ip, 1);
		df_dz.get(ip, z_size - 1) = df_dz.get(ip, z_size - 2);
	}
}

//Other prescriptions for derivative
//dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
//dfdz.get(ip, iz) = fabs(-fcr.get(ip, iz + 2) + 4.0 * fcr.get(ip, iz + 1) - 3.0 * fcr.get(ip, iz)) / 2. / dz;
//dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
