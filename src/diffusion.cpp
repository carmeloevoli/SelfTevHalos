// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "waves.h"

namespace CRWAVES {

void Waves::compute_D_zz() {
#pragma omp parallel for
  for (size_t ip = 0; ip < p_size; ++ip) {
    double k = 1. / larmor_radius(p.at(ip), par.magnetic_field);
    double factor = beta(p.at(ip)) * mks::c_light / 3.0 / pow2(k);
    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      const double value = factor / W_sg.get(ip, iz);
      D_zz.get(ip, iz) = value;
    }
    D_zz.get(ip, 0) = D_zz.get(ip, 1);
    D_zz.get(ip, z_size - 1) = D_zz.get(ip, z_size - 2);
  }
}

void Waves::compute_dfdz() {
  if (par.do_3D)
    compute_dfdz_3D();
  else
    compute_dfdz_1D();
}

void Waves::compute_dfdz_1D() {
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

void Waves::compute_dfdz_3D() {
  for (size_t ip = 0; ip < p_size; ++ip) {
    double value = 0;
    for (size_t iz = 0; iz < z_size - 1; ++iz) {
      if (iz == 0)
        value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
      else
        // value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
        value = fabs((f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1)) / 2. / dz);
      df_dz.get(ip, iz) = fabs(value);
    }
    df_dz.get(ip, 0) = df_dz.get(ip, 1);
    df_dz.get(ip, z_size - 1) = df_dz.get(ip, z_size - 2);
  }
}

// Other prescriptions for derivative tested
// dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
// dfdz.get(ip, iz) = fabs(-fcr.get(ip, iz + 2) + 4.0 * fcr.get(ip, iz + 1) - 3.0 * fcr.get(ip, iz)) / 2. / dz;
// dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;

}  // namespace CRWAVES