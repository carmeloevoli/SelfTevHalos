// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef PARAMS_H_
#define PARAMS_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "common.h"
#include "units.h"

namespace CRWAVES {

class Params {
 public:
  Params();
  virtual ~Params();
  void print();
  std::string generate_output_filename();

  void set_init_filename(std::string init_filename) { m_init_filename = init_filename; }
  void set_do_selfgeneration(bool do_selfgeneration) { m_do_selfgeneration = do_selfgeneration; }
  void set_alpha(double alpha) { m_alpha = alpha; };
  void set_magnetic_field(double magnetic_field) { m_magnetic_field = magnetic_field; };
  void set_do_kolmogorov(bool do_kolmogorov) { m_do_kolmogorov = do_kolmogorov; };
  void set_do_3D(bool do_3D) { m_do_3D = do_3D; };
  void set_spin_down_luminosity(double luminosity) { m_spin_down_luminosity = luminosity; }
  void set_source_size(double source_size) { m_source_size = source_size; }
  void set_p_size(size_t p_size) { m_p_size = p_size; }
  void set_z_size(size_t z_size) { m_z_size = z_size; }

  const size_t& p_size = m_p_size;
  const size_t& z_size = m_z_size;
  const double& correlation_length = m_correlation_length;
  const double& D_gal = m_D_gal;
  const double& D_gal_ref = m_D_gal_ref;
  const double& halo_size = m_halo_size;
  const double& p_min = m_p_min;
  const double& p_max = m_p_max;
  const double& magnetic_field = m_magnetic_field;
  const double& ck = m_ck;
  const double& alpha = m_alpha;
  const double& ion_number_density = m_ion_number_density;
  const double& source_pmin = m_source_pmin;
  const double& source_cutoff = m_source_cutoff;
  const double& source_size = m_source_size;
  const double& source_tdecay = m_source_tdecay;
  const double& spin_down_luminosity = m_spin_down_luminosity;
  const double& age = m_age;
  const double& tube_radius = m_tube_radius;
  const bool& do_selfgeneration = m_do_selfgeneration;
  const bool& do_3D = m_do_3D;
  const bool& do_kolmogorov = m_do_kolmogorov;
  const std::string& init_filename = m_init_filename;

 private:
  size_t m_p_size;
  size_t m_z_size;
  double m_correlation_length;
  double m_D_gal;
  double m_D_gal_ref;
  double m_halo_size;
  double m_p_min;
  double m_p_max;
  double m_magnetic_field;
  double m_ck;
  double m_alpha;
  double m_ion_number_density;
  double m_source_pmin;
  double m_source_cutoff;
  double m_source_size;
  double m_source_tdecay;
  double m_spin_down_luminosity;
  double m_age;
  double m_tube_radius;
  bool m_do_selfgeneration;
  bool m_do_3D;
  bool m_do_kolmogorov;
  std::string m_init_filename;
};

}  // namespace CRWAVES

#endif /* PARAMS_H_ */
