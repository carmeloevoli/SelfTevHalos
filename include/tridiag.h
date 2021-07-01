#include <cstdlib>
#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>

int solve_tridiag_nonsym(const std::vector<double> & diag, const std::vector<double> & abovediag, const std::vector<double> & belowdiag,
		const std::vector<double> & rhs, std::vector<double> & x, size_t N);

int gsl_linalg_solve_tridiag(const std::vector<double> & diag, const std::vector<double> & abovediag, const std::vector<double> & belowdiag,
		const std::vector<double> & rhs, std::vector<double> & solution);
