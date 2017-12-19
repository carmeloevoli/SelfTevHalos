#ifndef PARAMS_H_
#define PARAMS_H_

#include <iostream>
#include "units.h"
#include "utilities.h"

class Params {
private:
	double _correlation_length;
	double _k0;
	double _D_gal;
	double _D_gal_ref;
	double _vA_infty;
	double _magnetic_field;
	double _magnetic_energy_density;
	double _ck;
	double _alpha;
	double _ion_number_density;
	double _source_pmin;
	double _source_cutoff;
	double _source_size;
	double _source_tdecay;
	double _spin_down_luminosity;
	double _age;
	bool _do_selfgeneration;
	bool _do_3d;
	std::string _init_filename;

public:
	Params();
	std::string init_filename() const {
		return _init_filename;
	}
	double alpha() const {
		return _alpha;
	}
	double ck() const {
		return _ck;
	}
	double magnetic_field() const {
		return _magnetic_field;
	}
	double magnetic_energy_density() const {
		return _magnetic_energy_density;
	}
	double k0() const {
		return _k0;
	}
	double correlation_length() const {
		return _correlation_length;
	}
	double source_cutoff() const {
		return _source_cutoff;
	}
	double source_pmin() const {
		return _source_pmin;
	}
	double source_size() const {
		return _source_size;
	}
	double source_tdecay() const {
		return _source_tdecay;
	}
	double vA_infty() const {
		return _vA_infty;
	}
	void set_vA_infty(const double& vA) {
		_vA_infty = vA;
	}
	double spin_down_luminosity() const {
		return _spin_down_luminosity;
	}
	double age() const {
		return _age;
	}
	double D_gal() const {
		return _D_gal;
	}
	double D_gal_ref() const {
		return _D_gal_ref;
	}
	bool do_selfgeneration() const {
		return _do_selfgeneration;
	}
	bool do_3D() const {
		return _do_3d;
	}
	void set_selfgeneration() {
		_do_selfgeneration = true;
	}
	void unset_selfgeneration() {
		_do_selfgeneration = false;
	}
	void print();
};

#endif /* PARAMS_H_ */
