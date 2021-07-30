// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

#define sgn(A) ((A < 0.) ? -1. : 1.)

namespace CRWAVES {

#define pow2 utils::pow_integer<2>

Waves::Waves(const Params& params) : par(params) {
  magnetic_energy_density = pow2(par.magnetic_field) / 2. / mks::vacuum_permeability;
  vA_infty = par.magnetic_field / std::sqrt(mks::vacuum_permeability * mks::proton_mass * par.ion_number_density);
  k0 = 1. / par.correlation_length;
  factor_damping = par.ck * vA_infty;
  factor_growth = 2. * M_PI / 3. * mks::c_light * vA_infty / magnetic_energy_density;

  Q_cr.set_grid_size(par.p_size, par.z_size);
  W_ISM.set_grid_size(par.p_size, par.z_size);
  W_sg.set_grid_size(par.p_size, par.z_size);
  D_zz.set_grid_size(par.p_size, par.z_size);
  f_cr.set_grid_size(par.p_size, par.z_size);
  df_dz.set_grid_size(par.p_size, par.z_size);
  dp_dt.set_grid_size(par.p_size, 1);
  v_A.set_grid_size(1, par.z_size);
}

Waves::~Waves() { std::cout << "calling main destructor\n"; }

void Waves::build_p_axis() {
  p = utils::buildLogAxis<double>(par.p_min, par.p_max, par.p_size);
  this->p_size = p.size();
  utils::infoAxis<double>(p, "p", mks::GeV_c);
  dlnp = std::log(p.at(1) / p.at(0));
}

void Waves::build_z_axis() {
  if (par.do_3D)
    z = utils::buildLinAxis<double>(0, par.halo_size, par.z_size);
  else
    z = utils::buildLinAxis<double>(-par.halo_size, par.halo_size, par.z_size);
  this->z_size = z.size();
  utils::infoAxis(z, "z", mks::kpc);
  dz = fabs(z.at(1) - z.at(0));
}

void Waves::build_W_ISM() {
  const double k_norm = 1. / larmor_radius(par.D_gal_ref, par.magnetic_field);
  const double eta_B = mks::c_light / 2. / k_norm / par.D_gal * pow(k_norm / k0, 2. / 3.);
  std::cout << " - eta_B = " << eta_B << "\n";
  for (size_t ip = 0; ip < p.size(); ++ip) {
    for (size_t iz = 0; iz < z.size(); ++iz) {
      const double k = 1. / larmor_radius(p.at(ip), par.magnetic_field);
      double value = 2.0 * eta_B / 3.0 / k0;
      if (k > k0) value *= pow(k / k0, -5. / 3.);
      W_ISM.get(ip, iz) = value;
    }
  }
  W_ISM.show_grid("W_ISM", 1.);
}

void Waves::build_W_sg() {
  for (size_t ip = 0; ip < p.size(); ++ip) {
    for (size_t iz = 0; iz < z.size(); ++iz) {
      W_sg.get(ip, iz) = W_ISM.get(ip, iz);
    }
  }
  W_sg.show_grid("W_sg", 1.);
}

void Waves::build_f_cr() {
  for (size_t ip = 0; ip < p.size(); ++ip) {
    for (size_t iz = 0; iz < z.size(); ++iz) {
      df_dz.get(ip, iz) = 0;
      f_cr.get(ip, iz) = 0;
    }
  }
  df_dz.show_grid("dfdz", 1.);
  f_cr.show_grid("fp", 1.);
}

void Waves::build_D_zz() {
  compute_D_zz();
  D_zz.show_grid("D_zz", 1.);
}

void Waves::build_v_A() {
  for (size_t j = 0; j < z.size(); ++j) {
    // v_A.get(j) = sgn(z.at(j)) * par.vA_infty();
    v_A.get(j) = std::tanh(z.at(j) / mks::pc) * vA_infty;
  }
  v_A.show_grid("v_A", mks::cm / mks::s);
}

// void Waves::build_analytical_solution() { solution.set_params(par, compute_constant_CR_source_term()); }

}  // namespace CRWAVES