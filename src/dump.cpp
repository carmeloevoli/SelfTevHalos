// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include <iomanip>

#include "waves.h"

namespace CRWAVES {

#define pow2 utils::pow_integer<2>
#define pow3 utils::pow_integer<3>

std::string Waves::generate_output_filename(const std::string& s, const double& t) {
  std::stringstream sstream;
  sstream << "output/" << s << "_" << par.init_filename;
  sstream << "_nz_" << par.z_size << "_np_" << par.p_size;
  sstream << "_t_" << t / cgs::kyr;
  sstream << ".txt";
  std::string out = sstream.str();
  return out;
}

void Waves::dump(const double& t) {
  std::string filename = generate_output_filename("fcr", t);
  std::cout << "dumping f_cr on this file: " << filename << " ... ";
  std::ofstream outfile(filename.c_str());
  outfile << "#z[pc] p[GV] f[GV^-3 cm^-3] Dzz[cm^2/s] dfdz[GV^-3 cm^-3 kpc^-1] W_sg[pc] Q_cr[mks] tau_l [s]\n";
  outfile << std::scientific << std::setprecision(5);
  for (size_t iz = 0; iz < z.size(); ++iz)
    for (size_t ip = 0; ip < p.size(); ++ip) {
      outfile << z.at(iz) / cgs::parsec << " ";
      outfile << p.at(ip) / cgs::GV << " ";
      outfile << f_cr.get(ip, iz) / (1. / pow3(cgs::GV) / pow3(cgs::cm)) << " ";
      outfile << D_zz.get(ip, iz) / (pow2(cgs::cm) / cgs::second) << " ";
      outfile << df_dz.get(ip, iz) / (1. / pow3(cgs::GV) / pow3(cgs::cm) / cgs::kpc) << " ";
      outfile << W_sg.get(ip, iz) / (cgs::parsec) << " ";
      outfile << Q_cr.get(ip, iz) / (1. / pow3(cgs::GV) / pow3(cgs::cm) / cgs::second) << " ";
      outfile << p.at(ip) / -dp_dt.at(ip) / cgs::kyr << " ";
      outfile << "\n";
    }
  outfile.close();
  std::cout << "... done!\n";
}

// void Waves::dump_single(const double& t, const double& z_dump, const double& p_dump) {
//   std::string filename = generate_output_filename("fcr", t);
//   std::cout << "dumping f_cr on this file: " << filename << " ... ";
//   std::ofstream outfile(filename.c_str());
//   outfile << "# z[pc]  p[GeV/c]  fcr[GeV^-3 m^-3]  Dzz[cm^2/s]  dfdz[GeV^-3 m^-3 kpc^-1]  diff.flux[GeV^-3 "
//              "m^-3 s^-1] W_sg[mks] Q_cr[mks]\n";
//   outfile << std::scientific << std::setprecision(4);
//   const auto iz = utils::closestIndex(z_dump, z);
//   const auto ip = utils::closestIndex(p_dump, p);
//   {
//     outfile << z.at(iz) / mks::pc << " ";
//     outfile << p.at(ip) / mks::GeV_c << " ";
//     outfile << f_cr.get(ip, iz) / (1. / pow3(mks::GeV_c) / pow3(mks::meter)) << " ";
//     outfile << D_zz.get(ip, iz) / (pow2(mks::cm) / mks::s) << " ";
//     outfile << df_dz.get(ip, iz) / (1. / pow3(mks::GeV_c) / pow3(mks::meter) / mks::kpc) << " ";
//     outfile << D_zz.get(ip, iz) * df_dz.get(ip, iz) / (1. / pow3(mks::GeV_c) / pow2(mks::meter) / mks::s) << " ";
//     outfile << W_sg.get(ip, iz) << " ";
//     outfile << Q_cr.get(ip, iz) << " ";
//     outfile << "\n";
//   }
//   outfile.close();
//   std::cout << "... done!\n";
// }

// void Waves::dump_analytical_test(const double& t) {
//   std::string filename = generate_output_filename("test", t);
//   std::cout << "dumping analytical solution on this file: " << filename << " ... ";
//   std::ofstream outfile(filename.c_str());
//   outfile << "# z[pc]  p[GeV/c]  fcr[GeV^-3 m^-3]  fcr_a[GeV^-3 m^-3] \n";
//   outfile << std::scientific << std::setprecision(4);
//   double units = (1. / pow3(mks::GeV_c) / pow3(mks::meter));
//   for (size_t iz = 0; iz < z.size(); ++iz)
//     for (size_t ip = 0; ip < p.size(); ++ip) {
//       outfile << z.at(iz) / mks::pc << " ";
//       outfile << p.at(ip) / mks::GeV_c << " ";
//       outfile << f_cr.get(ip, iz) / units << " ";
//       outfile << solution.f(z.at(iz), t, p.at(ip)) / units << " ";
//       outfile << "\n";
//     }
//   outfile.close();
//   std::cout << "... done!\n";
// }

}  // namespace CRWAVES
