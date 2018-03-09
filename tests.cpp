//void Waves::build_CR_source_term_test() {
//	Q_cr.set_grid_size(p.get_size(), z.get_size());
//	/*double H = z.get_max();
//	 for (size_t ip = 0; ip < p.get_size(); ++ip) {
//	 for (size_t iz = 0; iz < z.get_size(); ++iz) {
//	 double r = max(fabs(z.at(iz)), 5. * pc);
//	 double D0 = pc * pc / year;
//	 double COS = cos(M_PI * r / 2. / H);
//	 double SIN = sin(M_PI * r / 2. / H);
//	 double A = pow2(M_PI) / 4. * (1. + pow2(r / H));
//	 double B = M_PI * r / H + M_PI * H / r * (1. + pow2(r / H));
//	 Q_cr.get(ip, iz) = D0 / pow2(H) * (A * COS + B * SIN);
//	 }
//	 }*/
//	for (size_t ip = 0; ip < p.get_size(); ++ip) {
//		for (size_t iz = 0; iz < z.get_size(); ++iz) {
//			double Q0 = 1.;
//			Q_cr.get(ip, iz) = Q0 * pow(p.at(ip) / electron_mass_c, -3.5);
//		}
//	}
//	Q_cr.show_grid("Q_CR", 1.);
//}
//
//void Waves::build_energy_losses_test() {
//	dpdt.set_grid_size(p.get_size(), z.get_size());
//	for (size_t iz = 0; iz < z.get_size(); ++iz) {
//		for (size_t ip = 0; ip < p.get_size(); ++ip) {
//			double gamma_e = p.at(ip) / electron_mass_c;
//			double normalization = 1.02e-16 * GeV_c / s;
//			dpdt.get(ip, iz) = -normalization * pow2(gamma_e);
//		}
//	}
//	dpdt.show_grid("dpdt", 1.);
//}
//
//void Waves::init_analytical_solution() {
//	double Q0 = compute_constant_CR_source_term();
//	solution.set_Q0(Q0);
//	solution.set_D0(params.D_gal());
//	solution.set_D0_ref(params.D_gal_ref());
//	solution.set_delta(1. / 3.);
//	solution.set_alpha(params.alpha());
//}
//
//void Waves::dump_fcr_test(const double& t) {
//	double Q0 = 1;
//	double b0 = 1.02e-16 * GeV_c / s;
//	double p0 = electron_mass_c;
//
//	string filename = generate_output_filename_p("fcr", t);
//	cout << "dumping fcr on this file: " << filename << " ... ";
//	ofstream outfile(filename.c_str());
//	outfile << "# z [kpc] \t p [GeV/c] \t fcr [] \n";
//	outfile << scientific << setprecision(5);
//	for (size_t iz = 0; iz < z.get_size(); ++iz) {
//		for (size_t ip = 0; ip < p.get_size(); ++ip) {
//			outfile << z.at(iz) / pc << "\t";
//			outfile << p.at(ip) / GeV_c << "\t";
//			outfile << fcr.get(ip, iz) << "\t";
//			outfile << D_zz.get(ip, iz) << "\t";
//			outfile << 0 << "\t";
//			outfile << Q_cr.get(ip, iz) << "\t";
//			outfile << 0 << "\t";
//			double x = p.at(ip) / p0;
//			double xmax = p.get_max() / p0;
//			outfile << 2. * Q0 * p0 / b0 * pow(x, -4.) * (pow(x, -0.5) - pow(xmax, -0.5)) << "\t";
//			outfile << cos(M_PI * z.at(iz) / 2. / z.get_max()) << "\t";
//			outfile << "\n";
//		}
//	}
//	outfile.close();
//	cout << "... done!" << "\n";
//}
//
//void Waves::dump(const double& t) {
//	double factor_damping = pow(2. * params.ck(), -1.5) * params.vA_infty();
//	double factor_growth = 2. * M_PI / 3. * c_light * params.vA_infty() / params.magnetic_energy_density();
//
//	string filename = generate_output_filename_p("fcr", t);
//	cout << "dumping fcr on this file: " << filename << " ... ";
//	ofstream outfile(filename.c_str());
//	outfile << "# z [kpc] \t p [GeV/c] \t fcr [] \n";
//	outfile << scientific << setprecision(5);
//	for (size_t iz = 0; iz < z.get_size(); ++iz) {
//		for (size_t ip = 0; ip < p.get_size(); ++ip) {
//			outfile << z.at(iz) / pc << "\t";
//			outfile << p.at(ip) / GeV_c << "\t";
//			outfile << fcr.get(ip, iz) / (1. / pow3(GeV_c) / pow3(meter)) << "\t";
//			outfile << D_zz.get(ip, iz) / (pow2(cm) / s) << "\t";
//			outfile << dfdz.get(ip, iz) / (1. / pow3(GeV_c) / pow3(meter) / kpc) << "\t";
//			outfile << D_zz.get(ip, iz) * dfdz.get(ip, iz) / (1. / pow3(GeV_c) / pow2(meter) / s) << "\t";
//			//outfile << solution.dfdz(z.at(iz), t, p.at(ip)) / (1. / pow3(GeV_c) / pow3(meter) / kpc) << "\t";
//			//outfile << (pow2(0.25 * kpc) / D_zz.get(ip, iz)) / kyr << "\t";
//			//outfile << fabs(p.at(ip) / dpdt.get(ip, iz)) / kyr << "\t";
//			//double k_ = 1. / larmor_radius(p.at(ip), params.magnetic_field());
//			//outfile << W_sg.get(ip, iz) / (factor_growth / k_ * pow4(p.at(ip)) * dfdz.get(ip, iz)) / year << "\t"; // GROWTH
//			//outfile << 1. / (factor_damping * pow(k_, 1.5) * sqrt(W_sg.get(ip, iz))) / year << "\t"; // DAMPING
//			////outfile << ncr.get(ip, iz) * (1. / GeV / pow3(meter)) << "\t";
//			//outfile << c_light * pow3(p.at(ip)) * solution.f(z.at(iz), t, p.at(ip)) << "\t"; // m^-2 s^-1
//			//outfile << c_light * pow3(p.at(ip)) * fcr.get(ip, iz) << "\t"; // m^-2 s^-1
//			//outfile << Q_cr.get(ip, iz) / (1./ pow3(GeV_c) / pow3(meter) / s) << "\t";
//			//outfile << 4.0 * M_PI * pow2(z.at(iz)) * fcr.get(ip, iz) / compute_constant_CR_source_term() / spectrum(p.at(ip), params.alpha(), params.source_pmin(), params.source_cutoff()) / (t * params.source_tdecay() / (t + params.source_tdecay())) * pc << "\t";
//			outfile << "\n";
//		}
//	}
//	outfile.close();
//	cout << "... done!" << "\n";
//}


