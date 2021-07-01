// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#ifndef PARAMS_H_
#define PARAMS_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "units.h"
#include "utilities.h"

namespace CRWAVES {

template <typename T>
class Param {
 public:
  Param() {}
  virtual ~Param() {}
  T get() const { return value; }
  void set(const T &value) { this->value = value; }

 protected:
  T value = T();
};

class Params {
 public:
  Params();
  virtual ~Params();
  void print();
  std::string generate_output_filename();

  Param<double> correlation_length;
  Param<double> D_gal;
  Param<double> D_gal_ref;
  Param<double> magnetic_field;
  Param<double> ck;
  Param<double> alpha;
  Param<double> ion_number_density;
  Param<double> source_pmin;
  Param<double> source_cutoff;
  Param<double> source_size;
  Param<double> source_tdecay;
  Param<double> spin_down_luminosity;
  Param<double> age;
  Param<double> tube_radius;
  Param<bool> do_selfgeneration;
  Param<bool> do_3D;
  Param<bool> do_kolmogorov;
  Param<std::string> init_filename;
};

}  // namespace CRWAVES

#endif /* PARAMS_H_ */
