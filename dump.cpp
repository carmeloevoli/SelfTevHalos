#include "waves.h"

std::string Waves::generate_output_filename(const std::string& s, const double& t) {
	std::stringstream sstream;
	sstream << "output/" << s << "_" << params.init_filename();
	sstream << "_t_" << t / kyr << "_nz_" << z.get_size() << "_np_" << p.get_size() << ".txt";
	std::string out = sstream.str();
	return out;
}

void Waves::dump(const double& t) {
	std::string filename = generate_output_filename("fcr", t);
	std::cout << "dumping f_cr on this file: " << filename << " ... ";
	std::ofstream outfile(filename.c_str());
	outfile << "# z [kpc] \t p [GeV/c] \t fcr [] \n";
	outfile << std::scientific << std::setprecision(4);
	for (size_t iz = 0; iz < z.get_size(); ++iz)
		for (size_t ip = 0; ip < p.get_size(); ++ip) {
			outfile << z.at(iz) / pc << "\t";
			outfile << p.at(ip) / GeV_c << "\t";
			outfile << f_cr.get(ip, iz) / (1. / pow3(GeV_c) / pow3(meter)) << "\t";
			outfile << D_zz.get(ip, iz) / (pow2(cm) / s) << "\t";
			outfile << df_dz.get(ip, iz) / (1. / pow3(GeV_c) / pow3(meter) / kpc) << "\t";
			outfile << D_zz.get(ip, iz) * df_dz.get(ip, iz) / (1. / pow3(GeV_c) / pow2(meter) / s) << "\t";
			outfile << "\n";
		}
	outfile.close();
	std::cout << "... done!" << "\n";
}

void Waves::dump_rates(const double& t) {
	std::string filename = generate_output_filename("rates", t);
	std::cout << "dumping rates on this file: " << filename << " ... ";
	std::ofstream outfile(filename.c_str());
}
