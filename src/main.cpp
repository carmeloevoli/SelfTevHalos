// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include <iostream>

#include "params.h"
#include "utilities/units.h"
#include "waves.h"

int main() {
  CRWAVES::Params par;
  par.set_init_filename("fiducial");
  par.print();

  CRWAVES::Waves W(par);
  W.build_z_axis();
  W.build_p_axis();
  W.build_v_A();
  W.build_energy_losses();
  W.build_CR_source_term();
  W.build_f_cr();
  W.build_Q_W();
  W.build_W();
  W.build_D_zz();

  W.dump(0);
  W.print_status(0, time(NULL));

  std::cout << " ... evolve\n";

  const double time_step = 0.1 * cgs::year;
  const int kyr_counter = 10 * 1000;
  const int max_counter = 1000 * kyr_counter;
  const int dump_counter = 1 * kyr_counter;

  W.evolve(time_step, max_counter, dump_counter);
  return 0;
}
