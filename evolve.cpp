#include <cmath>
#include <ctime>
#include <iostream>

#include "TGrid2D.h"
#include "units.h"
#include "waves.h"

void Waves::print_counter2time(const int& max_counter, const int& dump_counter) {
	cout << "running for " << max_counter * dt / kyr << " kyr with dt = " << dt / kyr << " kyr\n";
}

int Waves::get_difftime(time_t start) {
	time_t end = time(NULL);
	return difftime(end, start);
}

void Waves::print_status(const size_t& counter, const time_t& start) {
	cout << "dump at " << counter * dt / kyr << " kyr after " << get_difftime(start) << " secs\n";
	cout << " - f in " << fcr.min_value() << " ... " << fcr.max_value() << "\n";
	cout << " - dfdz in " << dfdz.min_value() << " ... " << dfdz.max_value() << "\n";
	cout << " - D in " << D_zz.min_value() / (cm2 / s) << " ... " << D_zz.max_value() / (cm2 / s) << "\n";
}

double Waves::compute_total_energy_in_fcr() {
	double value = 0;
	double dz = fabs(z.at(1) - z.at(0));
	double dlnp = std::log(p.at(1) / p.at(0));
	for (size_t iz = 0; iz < z.get_size(); ++iz) {
		double I_p = 0;
		for (size_t ip = 0; ip < p.get_size() - 1; ++ip) {
			I_p += pow3(p.at(ip)) * fcr.get(ip, iz) * (p.at(ip) * c_light);
		}
		//value += pow2(z.at(iz)) * I_p;
		value += I_p;
	}
	//value *= pow2(4. * M_PI) * dz * dlnp;
	double R = 2 * pc;
	value *= 4.0 * pow2(M_PI) * pow2(R) * dz * dlnp;
	return value;
}

double Waves::compute_source_luminosity() {
	double value = 0;
	double dz = fabs(z.at(1) - z.at(0));
	double dlnp = std::log(p.at(1) / p.at(0));
	for (size_t iz = 0; iz < z.get_size(); ++iz) {
		double I_p = 0;
		for (size_t ip = 0; ip < p.get_size() - 1; ++ip) {
			I_p += pow3(p.at(ip)) * Q_cr.get(ip, iz) * (p.at(ip) * c_light);
		}
		//value += pow2(z.at(iz)) * I_p;
		value += I_p;
	}
	//value *= pow2(4. * M_PI) * dz * dlnp;
	double R = 2 * pc;
	value *= 4.0 * pow2(M_PI) * pow2(R) * dz * dlnp;
	return value;
}

void Waves::compute_total_energy(const size_t& counter, const double& dt) {
	cout << "source term: " << compute_source_luminosity() / (erg / s) << "\n";
	double t = (double) counter * dt;
	double t0 = params.source_tdecay();
	injected_now = compute_source_luminosity() * (t * t0 / (t + t0));
	system_now = compute_total_energy_in_fcr();
	cout << "injected: " << (injected_now - injected_before) / erg << "\n";
	cout << "system: " << (system_now - system_before) / erg << "\n";
	cout << "ratio: " << (system_now - system_before) / (injected_now - injected_before) << "\n";
	injected_before = injected_now;
	system_before = system_now;
	cout << "\n";
}

void Waves::test_boundary_conditions() {
	double value = 0;
	for (size_t ip = 0; ip < p.get_size(); ++ip) {
		value += fabs(fcr.get(ip, z_size - 1));
		value += fabs(fcr.get(ip, 0));
	}
	for (size_t iz = 0; iz < z.get_size(); ++iz)
		value += fabs(fcr.get(p_size - 1, iz));
	cout << "boundary check -- value at border : " << value << "\n";
}

void Waves::evolve_f(const double& dt, const int& max_counter, const int& dump_counter) {
	set_dt(dt);
	print_counter2time(max_counter, dump_counter);
	time_t start = time(NULL);
	size_t counter = 0;
	while (counter < max_counter) {
		counter++;
		evolve_f_diffusive(1, counter * dt);
		//evolve_f_adiabatic(2, counter * dt);
		if (counter > 10 * 10000 && params.do_selfgeneration()) {
			compute_dfdz();
			//evolve_waves_advectice(1);
			evolve_waves_noadvectice(1);
			compute_D_zz();
			//	smooth_D_zz(3);
		}
		if (counter % dump_counter == 0) {
			print_status(counter, start);
			compute_total_energy(counter, dt);
			test_boundary_conditions();
			dump_fcr(counter * dt);
			//dump_fcr_test(counter * dt);
		} // if
	} // while
	print_status(counter, start);
}

void Waves::evolve() {
	//evolve_f(0.2 * year, 340 * 5000, 5000);
	evolve_f(0.1 * year, 340 * 10000, 10000);
}

