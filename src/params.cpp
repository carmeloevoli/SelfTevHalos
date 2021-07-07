// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "params.h"

namespace CRWAVES {

Params::Params() {
  m_init_filename = "fiducial";
  m_correlation_length = 10. * mks::parsec;
  m_source_cutoff = 100. * mks::TeV_c;
  m_source_pmin = 1 * mks::GeV_c;
  m_source_size = 0.1 * mks::pc;
  m_source_tdecay = 10 * mks::kyr;
  m_magnetic_field = 1 * mks::microgauss;
  m_ion_number_density = 1 / mks::cm3;
  m_ck = 5.2e-2;
  m_alpha = 3.5;
  m_age = 340 * mks::kyr;
  m_spin_down_luminosity = 3.8e34 * mks::erg / mks::s;
  m_tube_radius = 1 * mks::pc;
  m_D_gal = 5e28 * mks::cm2 / mks::s;
  m_D_gal_ref = 3 * mks::GeV_c;
  m_do_selfgeneration = false;
  m_do_3D = false;
  m_do_kolmogorov = true;
}

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
  double magnetic_energy_density = pow2(magnetic_field) / 2. / mks::vacuum_permeability / (mks::erg / mks::cm3);
  double vA_infty = magnetic_field / std::sqrt(mks::vacuum_permeability * mks::proton_mass * ion_number_density);
  std::string filename = generate_output_filename();
  std::cout << "dumping params on this file: " << filename << " ... ";
  std::ofstream parfile(filename.c_str());
  parfile << "l0 : " << correlation_length / mks::kpc << " kpc\n";
  parfile << "B0 : " << magnetic_field / mks::muG << " muG\n";
  parfile << "D_Gal : " << D_gal / (mks::cm2 / mks::s) << " cm2/s\n";
  parfile << "D_Gal_pref : " << D_gal_ref / mks::GeV_c << " GeV_c\n";
  parfile << "ck : " << ck << "\n";
  parfile << "injection slope : " << alpha << "\n";
  parfile << "injection cutoff : " << source_cutoff / mks::TeV_c << " TeV/c\n";
  parfile << "injection pmin : " << source_pmin / mks::GeV_c << " GeV/c\n";
  parfile << "source size  : " << source_size / mks::pc << " pc\n";
  parfile << "source t_decay : " << source_tdecay / mks::kyr << " kyr\n";
  parfile << "source age : " << age / mks::kyr << " kyr\n";
  parfile << "luminosity : " << spin_down_luminosity / (mks::erg / mks::s) << " erg/s\n";
  parfile << "tube radius : " << tube_radius / mks::pc << " pc\n";
  parfile << "do 3D? " << do_3D << "\n";
  parfile << "do self-generation? " << do_selfgeneration << "\n";
  parfile << "do Kolmogorov? " << do_kolmogorov << "\n";
  parfile << "\n";
  parfile << "magnetic_energy_density : " << magnetic_energy_density << " erg/cm3\n";
  parfile << "vA_infty : " << vA_infty / (mks::km / mks::s) << " km/s\n";
  parfile << "k0 : " << 1. / correlation_length * mks::parsec << " pc-1\n";
  parfile.close();
  std::cout << "done!\n";
}

}  // namespace CRWAVES