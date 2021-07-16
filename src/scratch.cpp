void evolve_waves_in_z(const size_t& number_of_operators);
void evolve_f_in_z_explicit(const size_t& number_of_operators, const double& t_now);
void evolve_waves_in_z_1D(const size_t& number_of_operators);
void evolve_waves_in_z_3D(const size_t& number_of_operators);

void Waves::evolve_waves_in_z_3D(const size_t& number_of_operators) {
  double vA_abs_dz = abs(vA_infty) / abs(z.at(1) - z.at(0));

#pragma omp parallel for
  for (size_t ip = 0; ip < p_size - 1; ++ip) {
    std::vector<double> rhs(z_size - 1);
    std::vector<double> central_diagonal(z_size - 1);
    std::vector<double> upper_diagonal(z_size - 2);
    std::vector<double> lower_diagonal(z_size - 2);
    std::vector<double> W_up(z_size - 1);

    double k_ = 1. / larmor_radius(p.at(ip), par.magnetic_field);

    for (size_t iz = 0; iz < z_size - 1; ++iz) {
      double L, C, U;
      L = vA_abs_dz;
      C = vA_abs_dz;
      U = 0.;

      central_diagonal.at(iz) = 1. + dt_half * C;
      if (iz != 0) {
        lower_diagonal.at(iz - 1) = -dt_half * L;
      }
      if (iz != z_size - 2) {
        upper_diagonal.at(iz) = -dt_half * U;
      }
      if (iz == 0) {
        upper_diagonal.at(iz) += -dt_half * L;
      }

      rhs.at(iz) = W_sg.get(ip, iz) * (2. - central_diagonal.at(iz));
      rhs.at(iz) += (iz != 0) ? dt_half * L * W_sg.get(ip, iz - 1) : 0.;
      rhs.at(iz) += (iz != z_size - 2) ? dt_half * U * W_sg.get(ip, iz + 1) : 0;

      double Gamma_D = Gamma_Damping(k_, W_sg.get(ip, iz));
      double Gamma_D_gal = factor_damping * pow(k_, 1.5) * sqrt(W_ISM.get(ip, iz));
      double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * df_dz.get(ip, iz);
      double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ISM.get(ip, iz);
      rhs.at(iz) += dt * Q_w / (double)number_of_operators;
    }

    GSL::gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, W_up);

    for (size_t iz = 0; iz < z_size - 1; ++iz) {
      double value = W_up.at(iz);  // min(W_up.at(iz), 1. / k_);
      W_sg.get(ip, iz) = std::max(value, 0.);
    }
  }
}

void Waves::evolve_waves_in_z(const size_t& number_of_operators) {
  if (par.do_3D)
    evolve_waves_in_z_3D(number_of_operators);
  else
    evolve_waves_in_z_1D(number_of_operators);
}

void Waves::evolve_waves_in_z_1D(const size_t& number_of_operators) {
  double vA_abs_dz = abs(vA_infty) / abs(z.at(1) - z.at(0));

#pragma omp parallel for
  for (size_t ip = 0; ip < p_size - 1; ++ip) {
    std::vector<double> rhs(z_size - 2);
    std::vector<double> central_diagonal(z_size - 2);
    std::vector<double> upper_diagonal(z_size - 3);
    std::vector<double> lower_diagonal(z_size - 3);
    std::vector<double> W_up(z_size - 2);

    double k_ = 1. / larmor_radius(p.at(ip), par.magnetic_field);

    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      double L, C, U;
      if (iz > z.get_idx()) {
        L = vA_abs_dz;
        C = vA_abs_dz;
        U = 0;
      } else if (iz < z.get_idx()) {
        L = 0;
        C = vA_abs_dz;
        U = vA_abs_dz;
      } else {
        L = -vA_abs_dz / 2.;
        C = 0;
        U = -vA_abs_dz / 2.;
      }

      central_diagonal.at(iz - 1) = 1. + dt_half * C;
      if (iz != 1) {
        lower_diagonal.at(iz - 2) = -dt_half * L;
      }
      if (iz != z_size - 2) {
        upper_diagonal.at(iz - 1) = -dt_half * U;
      }

      rhs.at(iz - 1) = W_sg.get(ip, iz) * (1. - dt_half * C);
      rhs.at(iz - 1) += (iz != 0) ? dt_half * L * W_sg.get(ip, iz - 1) : 0.;
      rhs.at(iz - 1) += (iz != z_size - 2) ? dt_half * U * W_sg.get(ip, iz + 1) : 0;

      double Gamma_D = Gamma_Damping(k_, W_sg.get(ip, iz));
      double Gamma_D_gal = factor_damping * pow(k_, 1.5) * sqrt(W_ISM.get(ip, iz));
      double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * df_dz.get(ip, iz);
      double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ISM.get(ip, iz);
      rhs.at(iz - 1) += dt * Q_w / (double)number_of_operators;
    }

    GSL::gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, W_up);

    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      double value = W_up.at(iz - 1);  // min(W_up.at(iz), 1. / k_);
      W_sg.get(ip, iz) = std::max(value, 0.);
    }
  }
}

void Waves::evolve_f_in_z_explicit(const size_t& number_of_operators, const double& t_now) {
  double vA_abs_dz = abs(vA_infty) / abs(z.at(1) - z.at(0));
  double dt_dz2 = dt / pow2(dz);

#pragma omp parallel for
  for (int ip = 0; ip < p_size - 1; ++ip) {
    std::vector<double> fcr_up(z_size, 0.0);
    for (int iz = 1; iz < z_size - 1; ++iz) {
      double Dz = D_zz.get(ip, iz);
      double DzUp = (iz < z_size - 1) ? D_zz.get(ip, iz + 1) : Dz;
      double DzDo = (iz > 0) ? D_zz.get(ip, iz - 1) : Dz;

      double UZ = 0.5 + (0.25 * (DzUp - DzDo) + Dz) * dt_dz2;
      double CZ = 2.0 * Dz * dt_dz2;
      double LZ = 0.5 - (0.25 * (DzUp - DzDo) - Dz) * dt_dz2;
      double Q = Q_cr.get(ip, iz) * source_evolution(t_now, par.source_tdecay);

      // fcr_up.at(iz) = UZ * f_cr.get(ip, iz + 1) - CZ * f_cr.get(ip, iz) + LZ * f_cr.get(ip, iz - 1);
      // fcr_up.at(iz) += dt * Q / (double) number_of_operators;
      fcr_up.at(iz) = f_cr.get(ip, iz);
      fcr_up.at(iz) += dt * Q / (double)number_of_operators;
      fcr_up.at(iz) += 0.25 * dt_dz2 * (DzUp - DzDo) * (f_cr.get(ip, iz + 1) - f_cr.get(ip, iz - 1));
      fcr_up.at(iz) += dt_dz2 * Dz * (f_cr.get(ip, iz + 1) - 2.0 * f_cr.get(ip, iz) + f_cr.get(ip, iz - 1));
    }
    for (int iz = 1; iz < z_size - 1; ++iz) {
      f_cr.get(ip, iz) = std::max(fcr_up.at(iz), 0.);
    }
  }
}

void Waves::evolve_waves_in_z(const size_t& number_of_operators) {
  double vA_abs_dz = abs(vA_infty) / abs(z.at(1) - z.at(0));

#pragma omp parallel for
  for (size_t ip = 0; ip < p_size - 1; ++ip) {
    std::vector<double> rhs(z_size - 2);
    std::vector<double> central_diagonal(z_size - 2);
    std::vector<double> upper_diagonal(z_size - 3);
    std::vector<double> lower_diagonal(z_size - 3);
    std::vector<double> W_up(z_size - 2);

    double k_ = 1. / larmor_radius(p.at(ip), par.magnetic_field);

    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      double L, C, U;
      if (iz > z.get_idx()) {
        L = vA_abs_dz;
        C = vA_abs_dz;
        U = 0;
      } else if (iz < z.get_idx()) {
        L = 0;
        C = vA_abs_dz;
        U = vA_abs_dz;
      } else {
        L = -vA_abs_dz / 2.;
        C = 0;
        U = -vA_abs_dz / 2.;
      }

      central_diagonal.at(iz - 1) = 1. + dt_half * C;
      if (iz != 1) {
        lower_diagonal.at(iz - 2) = -dt_half * L;
      }
      if (iz != z_size - 2) {
        upper_diagonal.at(iz - 1) = -dt_half * U;
      }

      rhs.at(iz - 1) = W_sg.get(ip, iz) * (1. - dt_half * C);
      rhs.at(iz - 1) += (iz != 0) ? dt_half * L * W_sg.get(ip, iz - 1) : 0.;
      rhs.at(iz - 1) += (iz != z_size - 2) ? dt_half * U * W_sg.get(ip, iz + 1) : 0;

      auto W = W_sg.get(ip, iz);
      auto Gamma_D = factor_damping * (par.do_kolmogorov) ? std::pow(k_, 1.5) * std::sqrt(W) : std::pow(k_, 2.) * W;
      auto Gamma_D_gal = factor_damping * std::pow(k_, 1.5) * std::sqrt(W_ISM.get(ip, iz));
      double WGamma_CR = factor_growth / k_ * pow4(p.at(ip)) * df_dz.get(ip, iz);
      double Q_w = WGamma_CR - Gamma_D * W_sg.get(ip, iz) + Gamma_D_gal * W_ISM.get(ip, iz);
      rhs.at(iz - 1) += dt * Q_w / (double)number_of_operators;
    }

    GSL::gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, W_up);

    for (size_t iz = 1; iz < z_size - 1; ++iz) {
      double value = W_up.at(iz - 1);  // min(W_up.at(iz), 1. / k_);
      W_sg.get(ip, iz) = std::max(value, 0.);
    }
  }
}