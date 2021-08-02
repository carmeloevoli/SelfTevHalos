// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include <cmath>
#include <ctime>
#include <iostream>

#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow3 utils::pow_integer<3>

int Waves::get_difftime(time_t start) const {
  time_t end = time(NULL);
  return difftime(end, start);
}

void Waves::print_status(const size_t& counter, const time_t& start) const {
  std::cout << "code status at " << (counter * dt) / cgs::kyr << " kyr\n";
  double f_units = 1. / pow3(cgs::GV) / pow3(cgs::cm);
  std::cout << " -> f in " << f_cr.min_value() / f_units << " ... " << f_cr.max_value() / f_units << "\n";
  std::cout << " -> dfdz in " << df_dz.min_value() << " ... " << df_dz.max_value() << "\n";
  double D_units = cgs::cm2 / cgs::second;
  std::cout << " -> D in " << D_zz.min_value() / D_units << " ... " << D_zz.max_value() / D_units << "\n";
  std::cout << "elapsed time : " << get_difftime(start) << " s\n";
}

// double Waves::compute_total_energy_in_fcr() {
//   const double R = 2 * mks::pc;
//   double value = 0;
//   for (size_t iz = 0; iz < z.size(); ++iz) {
//     double I_p = 0;
//     for (size_t ip = 0; ip < p.size() - 1; ++ip) {
//       I_p += pow3(p.at(ip)) * f_cr.get(ip, iz) * (p.at(ip) * mks::c_light);
//     }
//     if (par.do_3D)
//       value += pow2(z.at(iz)) * I_p;
//     else
//       value += I_p;
//   }
//   if (par.do_3D)
//     value *= pow2(4. * M_PI) * dz * dlnp;
//   else
//     value *= 4.0 * pow2(M_PI) * pow2(R) * dz * dlnp;

//   return value;
// }

// double Waves::compute_source_luminosity() {
//   const double R = 2 * mks::pc;
//   double value = 0;
//   for (size_t iz = 0; iz < z.size(); ++iz) {
//     double I_p = 0;
//     for (size_t ip = 0; ip < p.size() - 1; ++ip) {
//       I_p += pow3(p.at(ip)) * Q_cr.get(ip, iz) * (p.at(ip) * mks::c_light);
//     }
//     if (par.do_3D)
//       value += pow2(z.at(iz)) * I_p;
//     else
//       value += I_p;
//   }
//   if (par.do_3D)
//     value *= pow2(4. * M_PI) * dz * dlnp;
//   else
//     value *= 4.0 * pow2(M_PI) * pow2(R) * dz * dlnp;

//   return value;
// }

// void Waves::test_total_energy(const size_t& counter, const double& dt) {
//   std::cout << " -- source term: " << compute_source_luminosity() / (mks::erg / mks::s) << "\n";
//   double t = (double)counter * dt;
//   double t0 = par.source_tdecay;
//   double injected = compute_source_luminosity() * (t * t0 / (t + t0));
//   double in_crs = compute_total_energy_in_fcr();
//   std::cout << " -- injected: " << injected / mks::erg << "\n";
//   std::cout << " -- in CRs: " << in_crs / mks::erg << " / " << in_crs / injected << "\n";
// }

void Waves::test_boundary_conditions() {
  double value = 0;
  for (size_t ip = 0; ip < p_size; ++ip) value += fabs(f_cr.get(ip, z_size - 1));
  for (size_t iz = 0; iz < z.size(); ++iz) value += fabs(f_cr.get(p_size - 1, iz));
  std::cout << " -- total CRs at border : " << value << "\n";
}

void Waves::test_courant_conditions() {
  const double dt_dz2 = dt / pow2(z.at(1) - z.at(0));
  double beta_max = 0;
  for (size_t ip = 1; ip < p_size - 1; ++ip)
    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      double beta = D_zz.get(ip, iz) * dt_dz2;
      beta_max = std::max(beta, beta_max);
    }
  std::cout << " -- max CFL for diffusion " << beta_max << "\n";
}

void Waves::evolve(const double& dt, const int& max_counter, const int& dump_counter) {
  this->dt = dt;
  const time_t start = time(NULL);
  size_t counter = 0;
  while (counter < max_counter) {
    counter++;
    evolve_f_in_z(2, counter * dt);
    evolve_f_in_p(2, counter * dt);
    if ((double)counter * dt > 0.1 * cgs::kyr && par.do_selfgeneration) {
      compute_dfdz();
      compute_Q_W();
      // evolve_waves_in_z(1);
      evolve_waves();
      compute_D_zz();
    }
    if (counter % dump_counter == 0) {
      print_status(counter, start);
      //       test_total_energy(counter, dt);
      test_boundary_conditions();
      test_courant_conditions();
      // dump(counter * dt);
      dump_single(counter * dt, 10. * cgs::parsec, 10. * cgs::TV);  // TODO remove this at the end!
    }
  }
}

}  // namespace CRWAVES