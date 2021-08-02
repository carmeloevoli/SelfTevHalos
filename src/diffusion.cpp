// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>

void Waves::compute_D_zz() {
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double r_L = larmor_radius(p[ip], par.magnetic_field);
    const double factor = pow2(r_L) * cgs::c_light / 3.0;
    for (size_t iz = 1; iz < z.size() - 1; ++iz) {
      D_zz.get(ip, iz) = factor / W_sg.get(ip, iz);
    }
    D_zz.get(ip, 0) = D_zz.get(ip, 1);
    D_zz.get(ip, z.size() - 1) = D_zz.get(ip, z.size() - 2);
  }
}

// void Waves::compute_dfdz() {
//   if (par.do_3D)
//     compute_dfdz_3D();
//   else
//     compute_dfdz_1D();
// }

// void Waves::compute_dfdz_1D() {
//   for (size_t ip = 0; ip < p_size; ++ip) {
//     double value = 0;
//     for (size_t iz = 1; iz < z_size - 1; ++iz) {
//       if (iz == (z_size - 1) / 2)
//         value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
//       else
//         value = fabs((f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1)) / 2. / dz);
//       df_dz.get(ip, iz) = fabs(value);
//     }
//     df_dz.get(ip, 0) = df_dz.get(ip, 1);
//     df_dz.get(ip, z_size - 1) = df_dz.get(ip, z_size - 2);
//   }
// }

// void Waves::compute_dfdz_3D() {
//   for (size_t ip = 0; ip < p_size; ++ip) {
//     double value = 0;
//     for (size_t iz = 0; iz < z_size - 1; ++iz) {
//       if (iz == 0)
//         value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
//       else
//         // value = fabs(f_cr.get(ip, iz + 1) - f_cr.get(ip, iz)) / dz;
//         value = fabs((f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1)) / 2. / dz);
//       df_dz.get(ip, iz) = fabs(value);
//     }
//     df_dz.get(ip, 0) = df_dz.get(ip, 1);
//     df_dz.get(ip, z_size - 1) = df_dz.get(ip, z_size - 2);
//   }
// }

// Other prescriptions for derivative tested
// dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;
// dfdz.get(ip, iz) = fabs(-fcr.get(ip, iz + 2) + 4.0 * fcr.get(ip, iz + 1) - 3.0 * fcr.get(ip, iz)) / 2. / dz;
// dfdz.get(ip, iz) = fabs(fcr.get(ip, iz + 1) - fcr.get(ip, iz)) / dz;

}  // namespace CRWAVES