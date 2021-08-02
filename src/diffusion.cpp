// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow4 utils::pow_integer<4>

void Waves::compute_D_zz() {
  for (size_t ip = 0; ip < p_size; ++ip) {
    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      D_zz.get(ip, iz) = W2D(W_sg.get(ip, iz), p.at(ip), par.magnetic_field);
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
  const double dz = z.at(1) - z.at(0);
  for (size_t ip = 0; ip < p_size; ++ip) {
    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      double value = 0;
      if (iz == (z_size - 1) / 2)
        value = (f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
      else
        value = (f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1)) / 2. / dz;
      df_dz.get(ip, iz) = std::fabs(value);
    }
    df_dz.get(ip, 0) = df_dz.get(ip, 1);
    df_dz.get(ip, z_size - 1) = df_dz.get(ip, z_size - 2);
  }
}

void Waves::compute_dfdz_3D() {
  // const double dz = z.at(1) - z.at(0);
  // for (size_t ip = 0; ip < p_size; ++ip) {
  //   double value = 0;
  //   for (size_t iz = 0; iz < z_size - 1; ++iz) {
  //     if (iz == 0)
  //       value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
  //     else
  //       // value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
  //       value = fabs((f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1)) / 2. / dz);
  //     df_dz.get(ip, iz) = fabs(value);
  //   }
  //   df_dz.get(ip, 0) = df_dz.get(ip, 1);
  //   df_dz.get(ip, z_size - 1) = df_dz.get(ip, z_size - 2);
  // }
}

void Waves::compute_Q_W() {
  const double vA_infty = alfvenSpeed(par.ion_number_density, par.magnetic_field);
  const double U_B = magnetic_energy_density(par.magnetic_field);
  const double factor_damping = par.ck * vA_infty;
  const double factor_growth = 2. * M_PI / 3. * cgs::c_light * vA_infty / U_B;
#pragma omp parallel for
  for (size_t ip = 0; ip < p_size; ++ip) {
    const double k = 1. / larmor_radius(p.at(ip), par.magnetic_field);
    for (size_t iz = 0; iz < z_size; ++iz) {
      const double W = W_sg.get(ip, iz);
      const double W_ism = W_ISM.get(ip, iz);
      const double WGamma_CR = factor_growth / k * pow4(p.at(ip)) * df_dz.get(ip, iz);
      auto &Q = Q_W.get(ip, iz);
      Q = WGamma_CR;
      if (par.do_kolmogorov)
        Q += factor_damping * (-std::pow(k * W, 1.5) + std::pow(k * W_ism, 1.5));
      else {
        Q += factor_damping * (-std::pow(k * W, 2.) + std::pow(k * W_ism, 1.5));
      }
    }
  }
}

// Other prescriptions for derivative tested
// dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
// dfdz.get(ip, iz) = fabs(-fcr.get(ip, iz + 2) + 4.0 * fcr.get(ip, iz + 1) - 3.0 * fcr.get(ip, iz)) / 2. / dz;
// dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;

}  // namespace CRWAVES