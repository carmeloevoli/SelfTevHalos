#include "params.h"

Params::Params() {
	init_filename.set("pulsar_1D_5");
	correlation_length.set(10. * parsec);
	k0.set(1. / correlation_length.get());
	source_cutoff.set(100. * TeV_c);
	source_pmin.set(1 * GeV_c);
	source_size.set(0.1 * pc);
	source_tdecay.set(10 * kyr);
	magnetic_field.set(1.0 * microgauss);
	magnetic_energy_density.set(pow2(magnetic_field.get()) / 2. / vacuum_permeability);
	ion_number_density.set(1 / cm3);
	vA_infty.set(magnetic_field.get() / sqrt(vacuum_permeability * proton_mass * ion_number_density.get()));
	ck.set(5.2e-2);
	alpha.set(3.5);
	age.set(340 * kyr);
	spin_down_luminosity.set(3.8e34 * erg / s);
	D_gal.set(5e28 * cm2 / s);
	D_gal_ref.set(3 * GeV_c);
	tube_radius.set(1 * pc);
	do_selfgeneration.set(false);
	do_3D.set(false);
	do_kolmogorov.set(true);
}

Params::~Params() {
}

std::string Params::generate_output_filename() {
	std::stringstream sstream;
	sstream << "output/" << "par" << "_" << init_filename.get();
	sstream << ".txt";
	std::string out = sstream.str();
	return out;
}

void Params::print() {
	std::string filename = generate_output_filename();
	std::cout << "dumping params on this file: " << filename << " ... ";
	std::ofstream parfile(filename.c_str());
	//std::cout << "filename = " << _init_filename << "\n";
	parfile << "l0 = " << correlation_length.get() / kpc << " kpc \n";
	parfile << "k0 = " << k0.get() / (1. / kpc) << " kpc-1 \n";
	parfile << "vA = " << vA_infty.get() / (km / s) << " km/s \n";
	parfile << "B0 = " << magnetic_field.get() / muG << " muG \n";
	parfile << "U0 = " << magnetic_energy_density.get() / (eV / pow3(cm)) << " eV cm-3 \n";
	parfile << "D_Gal = " << D_gal.get() / (cm2 / s) << " cm2/s \n";
	parfile << "D_Gal_pref = " << D_gal_ref.get() / GeV_c << "\n";
	parfile << "ck = " << ck.get() << "\n";
	parfile << "injection slope = " << alpha.get() << "\n";
	parfile << "injection cutoff = " << source_cutoff.get() / TeV_c << " TeV/c \n";
	parfile << "injection pmin = " << source_pmin.get() / GeV_c << " GeV/c \n";
	parfile << "source size  = " << source_size.get() / pc << " pc \n";
	parfile << "source t_decay = " << source_tdecay.get() / kyr << " kyr \n";
	parfile << "source age = " << age.get() / kyr << " kyr \n";
	parfile << "luminosity = " << spin_down_luminosity.get() / (erg / s) << " erg/s \n";
	parfile << "tube radius = " << tube_radius.get() / pc << " pc\n";
	parfile << "do 3D? " << do_3D.get() << "\n";
	parfile << "do self-generation? " << do_selfgeneration.get() << "\n";
	parfile << "do Kolmogorov? " << do_kolmogorov.get() << "\n";
	parfile << "\n";
	parfile.close();
	std::cout << "done!\n";
}

