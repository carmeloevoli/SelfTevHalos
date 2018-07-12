#include <iostream>

#include "params.h"
#include "units.h"
#include "waves.h"

int main() {

	Params par;
	par.init_filename.set("new_geminga_7");
	par.do_selfgeneration.set(true);
	par.alpha.set(3.2);
	par.magnetic_field.set(2. * microgauss);
	par.do_kolmogorov.set(false);
	//par.do_3D.set(true);

	Waves * W = new Waves(par);

	std::cout << " ... init \n";

	W->build_z_axis(100 * pc, 1201);
	W->build_p_axis(1e2 * GeV_c, 1e6 * GeV_c, 64 * 4);
	W->build_CR_source_term();
	W->build_W_ISM();
	W->build_W_sg();
	W->build_f_cr();
	W->build_D_zz();
	W->build_v_A();
	W->build_energy_losses();
	//W->build_analytical_solution();

	W->dump(0);
	W->print_status(0, time(NULL));

	std::cout << " ... evolve\n";

	W->evolve(0.1 * year, 300 * 1000 * 10, 1000 * 10);

	std::cout << " ... done! \n";

	delete W;
	return 0;
}

