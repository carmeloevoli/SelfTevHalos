// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef WAVES_H_
#define WAVES_H_

// #include <algorithm>
// #include <cassert>
// #include <cmath>
// #include <cstdlib>
// #include <ctime>
// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <sstream>
// #include <vector>

#include "TGrid2D.h"
// #include "common.h"
#include "params.h"
#include "tridiag.h"
// #include "units.h"

namespace CRWAVES {

class Waves {
 public:
  Waves(const Params& params);
  virtual ~Waves();

  // builders.cpp
  void build_p_axis();
  void build_z_axis();
  void build_W();
  void build_f_cr();
  void build_Q_W();
  void build_D_zz();
  void build_v_A();
  void build_energy_losses();
  void build_CR_source_term();

  // dump.cpp
  std::string generate_output_filename(const std::string& s, const double& t);
  void dump(const double& t);
  void dump_single(const double& t, const double& z_dump, const double& p_dump);
  // // void dump_analytical_test(const double& t);

  // diffusion.cpp
  void compute_D_zz();
  void compute_dfdz();
  void compute_dfdz_1D();
  void compute_dfdz_3D();
  void compute_Q_W();

  // evolutors.cpp
  void evolve_f_in_z(const size_t& number_of_operators, const double& t_now);
  void evolve_f_in_z_1D(const size_t& number_of_operators, const double& t_now);
  void evolve_f_in_z_3D(const size_t& number_of_operators, const double& t_now);
  void evolve_f_in_p(const size_t& number_of_operators, const double& t_now);
  void evolve_waves();

  // // evolve.cpp
  // void set_dt(const double& dt);
  // void print_counter2time(const int& max_counter, const int& dump_counter);
  int get_difftime(time_t start) const;
  void print_status(const size_t& counter, const time_t& start) const;
  void compute_total_energy_in_fcr();
  void compute_source_luminosity(double t);
  // void test_total_energy(const size_t& counter, const double& dt);
  void test_boundary_conditions();
  void test_courant_conditions();
  void evolve(const double& dt, const int& max_counter, const int& dump_counter);

 private:
  Params par;

  size_t z_size = 0;
  size_t p_size = 0;
  double dt = 0;

  std::vector<double> p;
  std::vector<double> z;
  std::vector<double> v_A;
  std::vector<double> dp_dt;

  TGrid2D<double> Q_cr;
  TGrid2D<double> Q_W;
  TGrid2D<double> W_ISM;
  TGrid2D<double> W_sg;
  TGrid2D<double> f_cr;
  TGrid2D<double> df_dz;
  TGrid2D<double> D_zz;
};

}  // namespace CRWAVES

#endif /* WAVES_H_ */
