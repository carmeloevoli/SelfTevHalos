// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow3 utils::pow_integer<3>

void Waves::build_p_axis() {
  p = utils::buildLogAxis<double>(par.p_min, par.p_max, par.p_size);
  utils::infoAxis<double>(p, "p", cgs::GV);
  p_size = p.size();
}

void Waves::build_z_axis() {
  if (par.do_3D)
    z = utils::buildLinAxis<double>(0, par.halo_size, par.z_size);
  else
    z = utils::buildLinAxis<double>(-par.halo_size, par.halo_size, par.z_size);
  utils::infoAxis(z, "z", cgs::parsec);
  z_size = z.size();
}

void Waves::build_v_A() {
  v_A.resize(z_size);
  const double vA_infty = alfvenSpeed(par.ion_number_density, par.magnetic_field);
  for (size_t i = 0; i < z.size(); ++i) {
    v_A.at(i) = std::tanh(z.at(i) / cgs::parsec) * vA_infty;
  }
  utils::infoAxis(v_A, "v_A", cgs::km / cgs::second);
}

void Waves::build_W() {
  W_ISM.set_grid_size(p_size, z_size);
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double D = par.D_ISM * std::pow(p[ip] / par.D_ISM_p0, 1. / 3.);
    const double W = D2W(D, p.at(ip), par.magnetic_field);
    for (size_t iz = 0; iz < z.size(); ++iz) {
      W_ISM.get(ip, iz) = W;
    }
  }
  W_ISM.show_grid("W_ISM", cgs::parsec);

  W_sg.set_grid_size(p_size, z_size);
  W_sg.copy(W_ISM.getGrid());
  W_sg.show_grid("W_sg", cgs::parsec);
}

void Waves::build_f_cr() {
  f_cr.set_grid_size(p_size, z_size);
  f_cr.fill(0.);
  f_cr.show_grid("fp", 1.);
  df_dz.set_grid_size(p_size, z_size);
  df_dz.fill(0.);
  df_dz.show_grid("dfdz", 1.);
}

void Waves::build_Q_W() {
  Q_W.set_grid_size(p_size, z_size);
  Q_W.fill(0.);
  Q_W.show_grid("fp", 1.);
}

void Waves::build_D_zz() {
  D_zz.set_grid_size(p_size, z_size);
  compute_D_zz();
  D_zz.show_grid("D_zz", cgs::cm2 / cgs::second);
}

void Waves::build_energy_losses() {
  dp_dt.resize(p_size);
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double gamma_e = p.at(ip) / cgs::me_c;
    dp_dt.at(ip) = momentum_loss_rate(gamma_e, par.magnetic_field);
  }
  utils::infoAxis(dp_dt, "dp_dt", cgs::GV / cgs::second);
}

void Waves::build_CR_source_term() {
  Q_cr.set_grid_size(p_size, z_size);
  const double L_0 = initial_luminosity(par.source_luminosity_today, par.source_age, par.source_tdecay);
  const double Q_0 = (par.do_3D)
                         ? source_term_norm_3D(par.p_min, par.source_cutoff, par.source_slope, L_0, par.tube_radius)
                         : source_term_norm_1D(par.p_min, par.source_cutoff, par.source_slope, L_0, par.tube_radius);
  for (size_t ip = 0; ip < p.size(); ++ip) {
    const double F = pl_with_cutoff(p.at(ip), par.source_slope, par.source_pmin, par.source_cutoff);
    for (size_t iz = 0; iz < z.size(); ++iz) {
      const double z_abs = fabs(z.at(iz));
      const double G = (par.do_3D) ? gaussian_3D(z_abs, par.source_size) : gaussian_1D(z_abs, par.source_size);
      Q_cr.get(ip, iz) = Q_0 * G * F;
    }
  }
  Q_cr.show_grid("Q_CR", 1. / pow3(cgs::GV) / pow3(cgs::cm) / cgs::second);
}

}  // namespace CRWAVES