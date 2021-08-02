// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define sgn(A) ((A < 0.) ? -1. : 1.)
#define pow2 utils::pow_integer<2>

Waves::Waves(const Params& params) : par(params) {
  // k0 = 1. / par.correlation_length;
  // factor_damping = par.ck * vA_infty;
  // factor_growth = 2. * M_PI / 3. * mks::c_light * vA_infty / magnetic_energy_density;
  std::cout << "... main class constructed\n";
}

Waves::~Waves() { std::cout << "... main class destructed\n"; }

void Waves::build_p_axis() {
  p = utils::buildLogAxis<double>(par.p_min, par.p_max, par.p_size);
  utils::infoAxis<double>(p, "p", cgs::GV);
  //  dlnp = std::log(p.at(1) / p.at(0));
}

void Waves::build_z_axis() {
  if (par.do_3D)
    z = utils::buildLinAxis<double>(0, par.halo_size, par.z_size);
  else
    z = utils::buildLinAxis<double>(-par.halo_size, par.halo_size, par.z_size);
  utils::infoAxis(z, "z", cgs::parsec);
  // dz = fabs(z.at(1) - z.at(0));
}

void Waves::build_v_A() {
  v_A.resize(par.z_size);
  const double vA_infty = alfvenSpeed(par.ion_number_density, par.magnetic_field);
  for (size_t i = 0; i < z.size(); ++i) {
    v_A.at(i) = std::tanh(z.at(i) / cgs::parsec) * vA_infty;
  }
  utils::infoAxis(v_A, "v_A", cgs::km / cgs::second);
}

void Waves::build_W() {
  W_ISM.set_grid_size(par.p_size, par.z_size);
  // const double k_norm = 1. / larmor_radius(par.D_gal_ref, par.magnetic_field);
  //   const double eta_B = mks::c_light / 2. / k_norm / par.D_gal * pow(k_norm / k0, 2. / 3.);
  //   std::cout << " - eta_B = " << eta_B << "\n";
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double D = par.D_ISM * std::pow(p[ip] / par.D_ISM_p0, 1. / 3.);
    const double r_L = larmor_radius(p[ip], par.magnetic_field);
    const double value = cgs::c_light * pow2(r_L) / 3. / D;
    for (size_t iz = 0; iz < z.size(); ++iz) {
      //       const double k = 1. / larmor_radius(p.at(ip), par.magnetic_field);
      //       double value = 2.0 * eta_B / 3.0 / k0;
      //       if (k > k0) value *= pow(k / k0, -5. / 3.);
      W_ISM.get(ip, iz) = value;
    }
  }
  W_ISM.show_grid("W_ISM", cgs::parsec);

  W_sg.set_grid_size(par.p_size, par.z_size);
  W_sg.copy(W_ISM.getGrid());
  W_sg.show_grid("W_sg", cgs::parsec);
}

void Waves::build_f_cr() {
  f_cr.set_grid_size(par.p_size, par.z_size);
  f_cr.fill(0.);
  f_cr.show_grid("fp", 1.);

  df_dz.set_grid_size(par.p_size, par.z_size);
  df_dz.fill(0.);
  df_dz.show_grid("dfdz", 1.);
}

void Waves::build_D_zz() {
  D_zz.set_grid_size(par.p_size, par.z_size);
  compute_D_zz();
  D_zz.show_grid("D_zz", cgs::cm2 / cgs::second);
}

// void Waves::build_analytical_solution() { solution.set_params(par, compute_constant_CR_source_term()); }

}  // namespace CRWAVES