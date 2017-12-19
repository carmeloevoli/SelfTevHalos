#include <iostream>

#include "params.h"
#include "units.h"
#include "waves.h"

int main() {
	Waves * W = new Waves();

	W->build_z_axis(100 * pc, 2001);
	W->build_p_axis(1e2 * GeV_c, 1e6 * GeV_c, 5 * 4); // 32);
	W->build_CR_source_term();
	W->build_W_ISM();
	W->build_W_sg();
	W->build_f_cr();
	W->build_D_zz();

	W->build_v_A();
	W->build_energy_losses();

	std::cout << "dump initial state ... \n";
	W->dump(0);
	W->print_status(0, time(NULL));

	//cout << "evolve ... \n";

	//W->evolve();

	//cout << "... done! \n";

	delete W;
	return 0;
}
