#ifndef WAVES_H_
#define WAVES_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "tridiag.h"
#include "utilities.h"
#include "TAxis.h"
#include "TGrid2D.h"
#include "units.h"
#include "analytical_solution.h"
#include "params.h"

class Waves {
public:
	Waves();
	virtual ~Waves();

	//waves.cpp
	void build_p_axis(const double& pc_min, const double& pc_max, const size_t& pc_size);
	void build_z_axis(const double& halo_size, const size_t& k_size);
	double compute_constant_CR_source_term_1D();
	double compute_constant_CR_source_term_3D();

	void build_CR_source_term();
	void build_W_ISM();
	void build_W_sg();
	void build_f_cr();
	void build_D_zz();
	void build_v_A();
	void build_energy_losses();

	//dump.cpp
	std::string generate_output_filename(const std::string& s, const double& t);
	void dump(const double& t);
	void dump_rates(const double& t);



	//	void compute_total_energy(const size_t& counter, const double& dt);
//	double compute_total_energy_in_fcr();
//	double compute_source_luminosity();
//	void init_analytical_solution();
//

	//	void build_CR_source_term_test();
//	//dump.cpp
//	std::string generate_output_filename_p(const std::string& s, const double& t);
//	std::string generate_output_filename_k(const std::string& s, const double& t);
//	void dump_source();
//	void dump_wave_source();
//	void dump_CR_source();
//	void dump_Dzz(const double& t);
//	void dump_fcr(const double& t);
//	void dump_va();
//	void dump_fcr_test(const double& t);
//
//	//evolve.cpp
//	void evolve(const double& dt, const int& max_counter, const int& dump_counter);
//    void evolve();
//    void evolve_full();
//    void evolve_f(const double& dt_min, const double& dt_max, const double& dt_factor, const size_t& dt_iterations);
//    void evolve_f(const double& dt, const int& max_counter, const int& dump_counter);
//    double source_evolution(const double& t_now);
//
//    void print_counter2time(const int& max_counter, const int& dump_counter);
//	int get_difftime(time_t start);
//	void print_status(const size_t& counter, const time_t& start);
//
//	//evolutors.cpp
//	void evolve_f_diffusive(const size_t& number_of_operators, const double& t_now);
//	void evolve_f_adiabatic(const size_t& number_of_operators, const double& t_now);
//	void evolve_waves_advectice(const size_t& number_of_operators);
//	void evolve_waves_noadvectice(const size_t& number_of_operators);
//
//	void evolve_f_adiabatic_2nd(const size_t& number_of_operators, const double& t_now);
//	void evolve_f_adiabatic_2(const size_t& number_of_operators, const double& t_now);
//
//	//diffusion.cpp
//	std::vector<double> interpolate_W_k(const size_t& iz);
//	void compute_D_zz();
//	void smooth_D_zz(const int& smoothing_radius);
//	void smooth_W(const int& smoothing_radius);
//
//	std::vector<double> interpolate_dfdz(const size_t& iz);
//	void compute_Gamma_CR();
//	void compute_WGamma_CR();
//	void test_boundary_conditions();
//	double compute_dfdz();
//	void compute_D_zz_test();
//
//	//getters and setters
//	inline double get_dt() const {
//		return dt;
//	}
//	inline void set_dt(const double& dt) {
//		this->dt = dt;
//	}

private:
	double dt = double();
	size_t z_size = size_t();
	size_t p_size = size_t();
	Params params;
//
//	double injected_now = 0;
//	double injected_before = 0;
//
//	double system_now = 0;
//	double system_before = 0;
//
	TAxis<double> p;
	TAxis<double> z;
	TGrid2D<double> Q_cr;
	TGrid2D<double> W_ISM;
	TGrid2D<double> W_sg;
	TGrid2D<double> f_cr;
	TGrid2D<double> df_dz;
//	TGrid2D<double> D_kk;
	TGrid2D<double> D_zz;
	TGrid2D<double> v_A;
	TGrid2D<double> dp_dt;
//	TGrid2D<double> Gamma_CR;
//	TGrid2D<double> WGamma_CR;
};

#endif /* WAVES_H_ */
