// Copyright (c) 2018 carmeloevoli distributed under the MIT License
#include <iostream>

#include "params.h"
#include "units.h"
#include "waves.h"

int main() {
  CRWAVES::Params par;
  par.set_init_filename("test_3.5_kol");
  par.set_do_selfgeneration(true);
  par.set_alpha(3.5);
  par.set_magnetic_field(1. * mks::microgauss);
  par.set_do_kolmogorov(true);
  par.set_do_3D(false);

  CRWAVES::Waves W(par);

  std::cout << " ... init \n";

  W.build_z_axis(100 * mks::pc, 601);
  W.build_p_axis(1e2 * mks::GeV_c, mks::PeV_c, 32 * 4);
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

  double time_step = 0.1 * mks::year;
  double max_counter = 10 * 1000 * 300;  // 300 kyr
  double dump_counter = 10 * 1000;       // 1 kyr;

  W.evolve(time_step, max_counter, dump_counter);

  std::cout << " ... done! \n";
  return 0;
}
