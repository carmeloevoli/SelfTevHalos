// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "params.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>

Params::Params() {}
Params::~Params() {}

std::string Params::generate_output_filename() {
  std::stringstream sstream;
  sstream << "output/"
          << "par"
          << "_" << init_filename;
  sstream << ".txt";
  std::string out = sstream.str();
  return out;
}

void Params::print() {
  std::string filename = generate_output_filename();
  std::cout << "dumping params on this file: " << filename << " ... ";
  std::ofstream parfile(filename.c_str());
  parfile << "Input Parameters\n";
  parfile << "[p size : " << p_size << "]\n";
  parfile << "[z size : " << z_size << "]\n";
  parfile << "[p min : " << p_min / cgs::GV << " GV]\n";
  parfile << "[p max : " << p_max / cgs::GV << " GV]\n";
  parfile << "[halo size : " << halo_size / cgs::parsec << " pc]\n";
  parfile << "[magnetic field : " << magnetic_field / cgs::muG << " muG]\n";
  parfile << " D_ISM : " << D_ISM / (cgs::cm2 / cgs::second) << " cm2/s]\n";
  parfile << " D_ISM_p0 : " << D_ISM_p0 / cgs::GV << " GV]\n";
  parfile << "[ck : " << ck << "]\n";
  parfile << "[injection slope : " << source_slope << "]\n";
  parfile << "[injection cutoff : " << source_cutoff / cgs::TV << " TV]\n";
  parfile << "[injection pmin : " << source_pmin / cgs::GV << " GV]\n";
  parfile << "[source t_decay : " << source_tdecay / cgs::kyr << " kyr]\n";
  parfile << "[source age : " << source_age / cgs::kyr << " kyr]\n";
  parfile << "[luminosity : " << source_luminosity_today / (cgs::erg / cgs::second) << " erg/s]\n";
  parfile << "[tube radius : " << tube_radius / cgs::parsec << " pc]\n";
  parfile << "[do 3D? " << std::boolalpha << do_3D << "]\n";
  parfile << "[do self-generation? " << std::boolalpha << do_selfgeneration << "]\n";
  parfile << "[do Kolmogorov? " << std::boolalpha << do_kolmogorov << "]\n";
  parfile << "Derived Parameters\n";
  double vA_infty = alfvenSpeed(ion_number_density, magnetic_field);
  parfile << "[v_A : " << vA_infty / (cgs::km / cgs::second) << " km/s]\n";
  double u_B = magnetic_energy_density(magnetic_field);
  parfile << "[magnetic_energy_density : " << u_B / (cgs::eV / cgs::cm3) << " eV/cm3]\n";
  double L_0 = initial_luminosity(source_luminosity_today, source_age, source_tdecay);
  parfile << "[initial luminosity : " << L_0 / (cgs::erg / cgs::second) << " erg/s]\n";
  parfile.close();
  std::cout << "done!\n";
}

}  // namespace CRWAVES