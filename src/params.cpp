// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "params.h"

namespace CRWAVES {

Params::Params() {
  init_filename.set("fiducial");
  correlation_length.set(10. * mks::parsec);
  source_cutoff.set(100. * mks::TeV_c);
  source_pmin.set(1 * mks::GeV_c);
  source_size.set(0.1 * mks::pc);
  source_tdecay.set(10 * mks::kyr);
  magnetic_field.set(1.0 * mks::microgauss);
  ion_number_density.set(1 / mks::cm3);
  ck.set(5.2e-2);
  alpha.set(3.5);
  age.set(340 * mks::kyr);
  spin_down_luminosity.set(3.8e34 * mks::erg / mks::s);
  tube_radius.set(1 * mks::pc);
  D_gal.set(5e28 * mks::cm2 / mks::s);
  D_gal_ref.set(3 * mks::GeV_c);
  do_selfgeneration.set(false);
  do_3D.set(false);
  do_kolmogorov.set(true);
}

Params::~Params() {}

std::string Params::generate_output_filename() {
  std::stringstream sstream;
  sstream << "output/"
          << "par"
          << "_" << init_filename.get();
  sstream << ".txt";
  std::string out = sstream.str();
  return out;
}

void Params::print() {
  std::string filename = generate_output_filename();
  std::cout << "dumping params on this file: " << filename << " ... ";
  std::ofstream parfile(filename.c_str());
  parfile << "l0 : " << correlation_length.get() / mks::kpc << " kpc \n";
  parfile << "B0 : " << magnetic_field.get() / mks::muG << " muG \n";
  parfile << "D_Gal : " << D_gal.get() / (mks::cm2 / mks::s) << " cm2/s \n";
  parfile << "D_Gal_pref : " << D_gal_ref.get() / mks::GeV_c << " GeV_c\n";
  parfile << "ck : " << ck.get() << "\n";
  parfile << "injection slope : " << alpha.get() << "\n";
  parfile << "injection cutoff : " << source_cutoff.get() / mks::TeV_c << " TeV/c \n";
  parfile << "injection pmin : " << source_pmin.get() / mks::GeV_c << " GeV/c \n";
  parfile << "source size  : " << source_size.get() / mks::pc << " pc \n";
  parfile << "source t_decay : " << source_tdecay.get() / mks::kyr << " kyr \n";
  parfile << "source age : " << age.get() / mks::kyr << " kyr \n";
  parfile << "luminosity : " << spin_down_luminosity.get() / (mks::erg / mks::s) << " erg/s \n";
  parfile << "tube radius : " << tube_radius.get() / mks::pc << " pc\n";
  parfile << "do 3D? " << do_3D.get() << "\n";
  parfile << "do self-generation? " << do_selfgeneration.get() << "\n";
  parfile << "do Kolmogorov? " << do_kolmogorov.get() << "\n";
  parfile << "\n";
  parfile << "magnetic_energy_density : "
          << pow2(magnetic_field.get()) / 2. / mks::vacuum_permeability / (mks::erg / mks::cm3) << " erg/cm3\n";
  parfile << "vA_infty : "
          << magnetic_field.get() / std::sqrt(mks::vacuum_permeability * mks::proton_mass * ion_number_density.get()) /
                 (mks::km / mks::s)
          << " km/s\n";
  parfile << "k0 : " << 1. / correlation_length.get() * mks::parsec << " pc-1\n";
  parfile.close();
  std::cout << "done!\n";
}

}  // namespace CRWAVES