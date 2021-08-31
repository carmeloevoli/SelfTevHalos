// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "omp.h"
#include "waves.h"

#define OMP_NUM_THREADS THREADS

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow4 utils::pow_integer<4>

void Waves::evolve_f_in_z(const size_t& number_of_operators, const double& t_now) {
  if (par.do_3D)
    evolve_f_in_z_3D(number_of_operators, t_now);
  else
    evolve_f_in_z_1D(number_of_operators, t_now);
}

void Waves::evolve_f_in_z_1D(const size_t& number_of_operators, const double& t_now) {
  const double dz = std::abs(z.at(1) - z.at(0));
  const double dt_half = 0.5 * dt;
  //  const double vA_over_dz = abs(vA_infty) / dz;

#pragma omp parallel for num_threads(OMP_NUM_THREADS)
  for (int ip = 0; ip < p_size - 1; ++ip) {
    std::vector<double> rhs(z_size - 2);
    std::vector<double> central_diagonal(z_size - 2);
    std::vector<double> upper_diagonal(z_size - 3);
    std::vector<double> lower_diagonal(z_size - 3);
    std::vector<double> fcr_up(z_size - 2);

    for (int iz = 1; iz < z_size - 1; ++iz) {
      const double Dz = D_zz.get(ip, iz);
      const double DzUp = (iz < z_size - 1) ? D_zz.get(ip, iz + 1) : Dz;
      const double DzDo = (iz > 0) ? D_zz.get(ip, iz - 1) : Dz;

      const double Dz_dz2 = Dz / pow2(dz);
      const double dDz_4_dz2 = (DzUp - DzDo) / 4. / pow2(dz);

      double UZ = Dz_dz2 + dDz_4_dz2;
      double CZ = 2. * Dz_dz2;
      double LZ = Dz_dz2 - dDz_4_dz2;

      // TODO add advection

      central_diagonal.at(iz - 1) = 1. + dt_half * CZ;
      if (iz != z_size - 2) {
        upper_diagonal.at(iz - 1) = -dt_half * UZ;
      }
      if (iz != 1) {
        lower_diagonal.at(iz - 2) = -dt_half * LZ;
      }

      rhs.at(iz - 1) = f_cr.get(ip, iz) * (2. - central_diagonal.at(iz - 1));
      if (iz != z_size - 2) {
        rhs.at(iz - 1) -= f_cr.get(ip, iz + 1) * upper_diagonal.at(iz - 1);
      }
      if (iz != 1) {
        rhs.at(iz - 1) -= f_cr.get(ip, iz - 1) * lower_diagonal.at(iz - 2);
      }

      const double Q_t = Q_cr.get(ip, iz) * source_evolution(t_now, par.source_tdecay);
      rhs.at(iz - 1) += dt * Q_t / (double)number_of_operators;
    }

    GSL::gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

    for (int iz = 1; iz < z_size - 1; ++iz) {
      const double value = 0.5 * (fcr_up.at(iz - 1) + fcr_up.at(z_size - 2 - iz));
      f_cr.get(ip, iz) = std::max(value, 0.);
    }
  }  // for
}

void Waves::evolve_f_in_z_3D(const size_t& number_of_operators, const double& t_now) {
  const double dz = std::abs(z.at(1) - z.at(0));
  const double dt_half = 0.5 * dt;

  //   double vA_abs_dz = abs(vA_infty) / abs(z.at(1) - z.at(0));

#pragma omp parallel for num_threads(OMP_NUM_THREADS)
  for (int ip = 0; ip < p_size - 1; ++ip) {
    std::vector<double> rhs(z_size - 1);
    std::vector<double> central_diagonal(z_size - 1);
    std::vector<double> upper_diagonal(z_size - 2);
    std::vector<double> lower_diagonal(z_size - 2);
    std::vector<double> fcr_up(z_size - 1);

    for (int iz = 0; iz < z_size - 1; ++iz) {
      const double Dz = D_zz.get(ip, iz);
      const double DzUp = (iz < z_size - 1) ? D_zz.get(ip, iz + 1) : Dz;
      const double DzDo = (iz > 0) ? D_zz.get(ip, iz - 1) : Dz;

      const double Dz_dz2 = Dz / pow2(dz);
      const double dDz_4_dz2 = (DzUp - DzDo) / 4. / pow2(dz);
      const double r = std::max(fabs(z.at(iz)), 0.1 * cgs::parsec);  // TODO smaller radius
      const double Dz_dz_r = Dz / dz / r;

      double UZ = Dz_dz2 + dDz_4_dz2 + Dz_dz_r;
      double CZ = 2. * Dz_dz2;
      double LZ = Dz_dz2 - dDz_4_dz2 - Dz_dz_r;

      central_diagonal.at(iz) = 1. + dt_half * CZ;
      if (iz != z_size - 2) {
        upper_diagonal.at(iz) = -dt_half * UZ;
      }
      if (iz == 0) {
        upper_diagonal.at(iz) += -dt_half * LZ;
      }
      if (iz != 0) {
        lower_diagonal.at(iz - 1) = -dt_half * LZ;
      }
      rhs.at(iz) = f_cr.get(ip, iz) * (2. - central_diagonal.at(iz));
      if (iz != z_size - 2) {
        rhs.at(iz) -= f_cr.get(ip, iz + 1) * upper_diagonal.at(iz);
      }
      if (iz != 0) {
        rhs.at(iz) -= f_cr.get(ip, iz - 1) * lower_diagonal.at(iz - 1);
      }
      double Q = Q_cr.get(ip, iz) * source_evolution(t_now, par.source_tdecay);
      rhs.at(iz) += dt * Q / (double)number_of_operators;
    }

    GSL::gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

    for (int iz = 0; iz < z_size - 1; ++iz) {
      double value = fcr_up.at(iz);
      f_cr.get(ip, iz) = std::max(value, 0.);
    }
  }  // for
}

void Waves::evolve_f_in_p(const size_t& number_of_operators, const double& t_now) {
  const double dt_half = 0.5 * dt;

#pragma omp parallel for num_threads(OMP_NUM_THREADS)
  for (int iz = 1; iz < z_size - 1; ++iz) {
    std::vector<double> rhs(p_size - 1);
    std::vector<double> central_diagonal(p_size - 1);
    std::vector<double> upper_diagonal(p_size - 2);
    std::vector<double> lower_diagonal(p_size - 2);
    std::vector<double> fcr_up(p_size - 1);

    // double dvadz = (v_A.get(iz + 1) - v_A.get(iz - 1)) / (z.at(iz + 1) - z.at(iz - 1));

    for (int ip = 0; ip < p_size - 1; ++ip) {
      const double pUp = p.at(ip + 1) - p.at(ip);

      double LP = 0;
      double CP = -dp_dt.at(ip) / pUp;
      double UP = -pow2(p.at(ip + 1) / p.at(ip)) * dp_dt.at(ip + 1) / pUp;

      // double b_i = 0;  // dvadz / 3.0 / (p.at(ip + 1) / p.at(ip) - 1.0);
      // double L = 0, C = 0, U = 0;

      // C += b_i;
      // U += b_i;

      // C -= dp_dt.at(ip) / (p.at(ip + 1) - p.at(ip));
      // U -= pow2(p.at(ip + 1) / p.at(ip)) * dp_dt.at(ip + 1) / (p.at(ip + 1) - p.at(ip));

      central_diagonal.at(ip) = 1. + dt_half * CP;
      if (ip != 0) {
        lower_diagonal.at(ip - 1) = -dt_half * LP;
      }
      if (ip != p_size - 2) {
        upper_diagonal.at(ip) = -dt_half * UP;
      }
      rhs.at(ip) = f_cr.get(ip, iz) * (2. - central_diagonal.at(ip));
      if (ip != p_size - 2) {
        rhs.at(ip) -= f_cr.get(ip + 1, iz) * upper_diagonal.at(ip);
      }
      if (ip != 0) {
        rhs.at(ip) -= f_cr.get(ip - 1, iz) * lower_diagonal.at(ip - 1);
      }
      const double Q_t = Q_cr.get(ip, iz) * source_evolution(t_now, par.source_tdecay);
      rhs.at(ip) += dt * Q_t / (double)number_of_operators;
    }

    GSL::gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

    for (int ip = 0; ip < p_size - 1; ++ip) {
      const double value = fcr_up.at(ip);
      f_cr.get(ip, iz) = std::max(value, 0.);
    }
  }  // for
}

void Waves::evolve_waves() {
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
  for (size_t ip = 0; ip < p_size; ++ip) {
    for (size_t iz = 0; iz < z_size; ++iz) {
      const double value = W_sg.get(ip, iz) + dt * Q_W.get(ip, iz);
      // value = std::min(value, 1. / k_); TODO check if Bohm limit is violated
      W_sg.get(ip, iz) = std::max(value, 0.);
    }
  }
}

}  // namespace CRWAVES