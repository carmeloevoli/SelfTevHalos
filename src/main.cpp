// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include <iostream>

#include "params.h"
#include "units.h"
#include "waves.h"

int main() {
  CRWAVES::Params par;
  par.set_init_filename("testN0.1pc_3.5_kol_3.8e33");
  par.set_do_selfgeneration(true);
  par.set_alpha(3.5);
  par.set_magnetic_field(1. * mks::microgauss);
  par.set_do_kolmogorov(true);
  par.set_do_3D(false);
  par.set_spin_down_luminosity(3.8e33 * mks::erg / mks::s);
  par.set_source_size(0.1 * mks::pc);
  par.set_p_size(32 * 4);
  par.set_z_size(401);

  CRWAVES::Waves W(par);

  std::cout << " ... init \n";

  W.build_z_axis();
  W.build_p_axis();
  W.build_CR_source_term();
  W.build_W_ISM();
  W.build_W_sg();
  W.build_f_cr();
  W.build_D_zz();
  W.build_v_A();
  W.build_energy_losses();
  // W->build_analytical_solution();

  W.dump(0);
  W.print_status(0, time(NULL));

  std::cout << " ... evolve\n";

  const double time_step = 0.1 * mks::year;
  const int kyr_counter = 10 * 1000;
  const int max_counter = 1000 * kyr_counter;
  const int dump_counter = 1 * kyr_counter;

  W.evolve(time_step, max_counter, dump_counter);

  std::cout << " ... done! \n";
  return 0;
}
