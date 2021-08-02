// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef PARAMS_H_
#define PARAMS_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "common.h"
#include "utilities/units.h"

namespace CRWAVES {

class Params {
 public:
  Params();
  virtual ~Params();
  void print();
  std::string generate_output_filename();

  void set_init_filename(std::string init_filename) {
    m_init_filename = init_filename;
  }
  void set_z_size(size_t z_size) {
    m_z_size = z_size;
  }
  void set_p_size(size_t p_size) {
    m_p_size = p_size;
  }
  void set_source_luminosity(double luminosity) {
    m_source_luminosity_today = luminosity;
  }
  void set_source_slope(double alpha) {
    m_source_slope = alpha;
  };
  void set_magnetic_field(double magnetic_field) {
    m_magnetic_field = magnetic_field;
  };
  void set_do_kolmogorov(bool do_kolmogorov) {
    m_do_kolmogorov = do_kolmogorov;
  };
  void set_do_3D(bool do_3D) {
    m_do_3D = do_3D;
  };
  void set_do_selfgeneration(bool do_selfgeneration) {
    m_do_selfgeneration = do_selfgeneration;
  }

  const size_t& z_size = m_z_size;
  const size_t& p_size = m_p_size;
  const double& halo_size = m_halo_size;
  const double& p_min = m_p_min;
  const double& p_max = m_p_max;
  const double& D_ISM = m_D_ISM;
  const double& D_ISM_p0 = m_D_ISM_p0;
  const double& ck = m_ck;
  const double& magnetic_field = m_magnetic_field;
  const double& ion_number_density = m_ion_number_density;
  const double& source_age = m_source_age;
  const double& source_slope = m_source_slope;
  const double& source_pmin = m_source_pmin;
  const double& source_cutoff = m_source_cutoff;
  const double& source_tdecay = m_source_tdecay;
  const double& source_luminosity_today = m_source_luminosity_today;
  const double& tube_radius = m_tube_radius;
  const bool& do_selfgeneration = m_do_selfgeneration;
  const bool& do_3D = m_do_3D;
  const bool& do_kolmogorov = m_do_kolmogorov;
  const std::string& init_filename = m_init_filename;

 private:
  size_t m_z_size{401};
  size_t m_p_size{32 * 4};
  double m_halo_size{100. * cgs::parsec};
  double m_p_min{1e2 * cgs::GV};
  double m_p_max{1. * cgs::PV};
  double m_source_age{340 * cgs::kyr};
  double m_source_slope{3.5};
  double m_source_pmin{1 * cgs::GV};
  double m_source_cutoff{100. * cgs::TV};
  double m_source_tdecay{10 * cgs::kyr};
  double m_source_luminosity_today{3.8e34 * cgs::erg / cgs::second};
  double m_tube_radius{1 * cgs::parsec};
  double m_D_ISM{5e28 * cgs::cm2 / cgs::second};
  double m_D_ISM_p0{3 * cgs::GV};
  double m_ck{5.2e-2};
  double m_magnetic_field{cgs::muG};
  double m_ion_number_density{1. / cgs::cm3};
  bool m_do_selfgeneration{true};
  bool m_do_3D{false};
  bool m_do_kolmogorov{true};
  std::string m_init_filename{"fiducial"};
};

}  // namespace CRWAVES

#endif /* PARAMS_H_ */
