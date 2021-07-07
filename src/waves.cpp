// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

#define sgn(A) ((A < 0.) ? -1. : 1.)

namespace CRWAVES {

Waves::Waves(const Params& params) : par(params) {
  par.print();

  magnetic_energy_density = pow2(par.magnetic_field) / 2. / mks::vacuum_permeability;
  vA_infty = par.magnetic_field / std::sqrt(mks::vacuum_permeability * mks::proton_mass * par.ion_number_density);
  k0 = 1. / par.correlation_length;
  factor_damping = par.ck * vA_infty;
  factor_growth = 2. * M_PI / 3. * mks::c_light * vA_infty / magnetic_energy_density;
}

Waves::~Waves() { std::cout << "calling main destructor\n"; }

void Waves::build_p_axis(const double& p_min, const double& p_max, const size_t& p_size) {
  p.build_log_axis(p_min, p_max, p_size);
  this->p_size = p.get_size();
  p.set_reference_value(sqrt(p_min * p_max));
  p.show_axis("p", mks::GeV_c);
  dlnp = std::log(p.at(1) / p.at(0));
}

void Waves::build_z_axis(const double& halo_size, const size_t& z_size) {
  if (par.do_3D)
    z.build_lin_axis(0, halo_size, z_size);
  else
    z.build_lin_axis(-halo_size, halo_size, z_size);
  this->z_size = z.get_size();
  z.set_reference_value(0);
  z.show_axis("z", mks::kpc);
  dz = fabs(z.at(1) - z.at(0));
}

void Waves::build_W_ISM() {
  W_ISM.set_grid_size(p.get_size(), z.get_size());
  double k_norm = 1. / larmor_radius(par.D_gal_ref, par.magnetic_field);
  double eta_B = mks::c_light / 2. / k_norm / par.D_gal * pow(k_norm / k0, 2. / 3.);
  std::cout << " - eta_B = " << eta_B << "\n";
  for (size_t ip = 0; ip < p.get_size(); ++ip) {
    for (size_t iz = 0; iz < z.get_size(); ++iz) {
      double value = 2.0 * eta_B / 3.0 / k0;
      double k = 1. / larmor_radius(p.at(ip), par.magnetic_field);
      if (k > k0) value *= pow(k / k0, -5. / 3.);
      W_ISM.get(ip, iz) = value;
    }
  }
  W_ISM.show_grid("W_ISM", 1.);
}

void Waves::build_W_sg() {
  W_sg.set_grid_size(p.get_size(), z.get_size());
  for (size_t ip = 0; ip < p.get_size(); ++ip) {
    for (size_t iz = 0; iz < z.get_size(); ++iz) {
      W_sg.get(ip, iz) = W_ISM.get(ip, iz);
    }
  }
  W_sg.show_grid("W_sg", 1.);
}

void Waves::build_f_cr() {
  df_dz.set_grid_size(p.get_size(), z.get_size());
  f_cr.set_grid_size(p.get_size(), z.get_size());
  for (size_t ip = 0; ip < p.get_size(); ++ip) {
    for (size_t iz = 0; iz < z.get_size(); ++iz) {
      df_dz.get(ip, iz) = 0;
      f_cr.get(ip, iz) = 0;
    }
  }
  df_dz.show_grid("dfdz", 1.);
  f_cr.show_grid("fp", 1.);
}

void Waves::build_D_zz() {
  D_zz.set_grid_size(p.get_size(), z.get_size());
  compute_D_zz();
  D_zz.show_grid("D_zz", 1.);
}

void Waves::build_v_A() {
  v_A.set_grid_size(1, z.get_size());
  for (size_t j = 0; j < z.get_size(); ++j) {
    // v_A.get(j) = sgn(z.at(j)) * par.vA_infty();
    v_A.get(j) = std::tanh(z.at(j) / mks::pc) * vA_infty;
  }
  v_A.show_grid("v_A", mks::km / mks::s);
}

void Waves::build_analytical_solution() { solution.set_params(par, compute_constant_CR_source_term()); }

}  // namespace CRWAVES