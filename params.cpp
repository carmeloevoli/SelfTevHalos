#include "params.h"

Params::Params() {
	_init_filename = "test_1D_0";
	_correlation_length = 10. * parsec;
	_k0 = 1. / _correlation_length;
	_source_cutoff = 100. * TeV_c;
	_source_pmin = 1 * GeV_c;
	_source_size = 1 * pc;
	_source_tdecay = 8 * kyr;
	_magnetic_field = 1.0 * microgauss;
	_magnetic_energy_density = pow2(_magnetic_field) / 2. / vacuum_permeability;
	_ion_number_density = 1. / cm3;
	_vA_infty = _magnetic_field / sqrt(vacuum_permeability * proton_mass * _ion_number_density);
	_ck = 5.2e-2;
	_alpha = 3.5;
	_age = 340 * kyr;
	_spin_down_luminosity = 3.8e34 * erg / s;
	_D_gal = 5e28 * cm2 / s;
	_D_gal_ref = 3 * GeV_c;
	_do_selfgeneration = false;
	_do_3d = false;
}

void Params::print() {
	std::cout << "l0 = " << _correlation_length / kpc << " kpc \n";
	std::cout << "k0 = " << _k0 / (1./kpc) << " kpc-1 \n";
	std::cout << "vA = " << _vA_infty / (km / s) << " km/s \n";
	std::cout << "B0 = " << _magnetic_field / muG << " muG \n";
	std::cout << "U0 = " << _magnetic_energy_density / (eV / pow3(cm)) << " eV cm-3 \n";
	std::cout << "D_Gal = " << _D_gal / (cm2 / s) << " cm2/s \n";
	std::cout << "D_Gal_pref = " << _D_gal_ref / GeV_c << "\n";
	std::cout << "ck = " << _ck << "\n";
	std::cout << "injection slope = " << _alpha << "\n";
	std::cout << "injection cutoff = " << _source_cutoff / TeV_c << " TeV/c \n";
	std::cout << "injection pmin = " << _source_pmin / GeV_c << " GeV/c \n";
	std::cout << "source size  = " << _source_size / pc << " pc \n";
	std::cout << "source t_decay = " << _source_tdecay / kyr << " kyr \n";
	std::cout << "source age = " << _age / kyr << " kyr \n";
	std::cout << "luminosity = " << _spin_down_luminosity / (erg / s) << " erg/s \n";
	std::cout << "\n";
}

