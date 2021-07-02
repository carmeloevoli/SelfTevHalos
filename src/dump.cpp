// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include "waves.h"

namespace CRWAVES {

std::string Waves::generate_output_filename(const std::string& s, const double& t) {
  std::stringstream sstream;
  sstream << "output/" << s << "_" << par.init_filename.get();
  sstream << "_t_" << t / mks::kyr << "_nz_" << z.get_size() << "_np_" << p.get_size() << ".txt";
  std::string out = sstream.str();
  return out;
}

void Waves::dump(const double& t) {
  std::string filename = generate_output_filename("fcr", t);
  std::cout << "dumping f_cr on this file: " << filename << " ... ";
  std::ofstream outfile(filename.c_str());
  outfile << "# z[pc]  p[GeV/c]  fcr[GeV^-3 m^-3]  Dzz[cm^2/s]  dfdz[GeV^-3 m^-3 kpc^-1]  diff.flux[GeV^-3 "
             "m^-3 s^-1] W_sg[mks] Q_cr[mks]\n";
  outfile << std::scientific << std::setprecision(4);
  for (size_t iz = 0; iz < z.get_size(); ++iz)
    for (size_t ip = 0; ip < p.get_size(); ++ip) {
      outfile << z.at(iz) / mks::pc << " ";
      outfile << p.at(ip) / mks::GeV_c << " ";
      outfile << f_cr.get(ip, iz) / (1. / pow3(mks::GeV_c) / pow3(mks::meter)) << " ";
      outfile << D_zz.get(ip, iz) / (pow2(mks::cm) / mks::s) << " ";
      outfile << df_dz.get(ip, iz) / (1. / pow3(mks::GeV_c) / pow3(mks::meter) / mks::kpc) << " ";
      outfile << D_zz.get(ip, iz) * df_dz.get(ip, iz) / (1. / pow3(mks::GeV_c) / pow2(mks::meter) / mks::s) << " ";
      outfile << W_sg.get(ip, iz) << " ";
      outfile << Q_cr.get(ip, iz) << " ";
      outfile << "\n";
    }
  outfile.close();
  std::cout << "... done!"
            << "\n";
}

void Waves::dump_analytical_test(const double& t) {
  std::string filename = generate_output_filename("test", t);
  std::cout << "dumping analytical solution on this file: " << filename << " ... ";
  std::ofstream outfile(filename.c_str());
  outfile << "# z[pc]  p[GeV/c]  fcr[GeV^-3 m^-3]  fcr_a[GeV^-3 m^-3] \n";
  outfile << std::scientific << std::setprecision(4);
  double units = (1. / pow3(mks::GeV_c) / pow3(mks::meter));
  for (size_t iz = 0; iz < z.get_size(); ++iz)
    for (size_t ip = 0; ip < p.get_size(); ++ip) {
      outfile << z.at(iz) / mks::pc << " ";
      outfile << p.at(ip) / mks::GeV_c << " ";
      outfile << f_cr.get(ip, iz) / units << " ";
      outfile << solution.f(z.at(iz), t, p.at(ip)) / units << " ";
      outfile << "\n";
    }
  outfile.close();
  std::cout << "... done!\n";
}

}  // namespace CRWAVES
