#include <iostream>
#include <functional>

#include "murlib/interpolation.h"
#include "murlib/approximation.h"
#include "murlib/matrix_solver.h"
#include "murlib/integration.h"
#include "murlib/miscellaneous.h"



void murlib::mnrp_for_polynoms(const std::function<double(double)> f, const int n, const double a, const double b, const int test_grid_size, const double* grid_xi, double* result_yi, double* error) {
	// МНРП!

	double* chebyshev_grid_xi = new double[n];
	double* chebyshev_grid_yi = new double[n];

	for (int i = 0; i < n; ++i) {
		chebyshev_grid_xi[i] = murlib::get_chebyshev_root(a, b, n, i);
		chebyshev_grid_yi[i] = f(chebyshev_grid_xi[i]);
	}

	double max_error = 0;
	for (int i = 0; i < test_grid_size; ++i) {
		result_yi[i] = murlib::lagrange_polynom(chebyshev_grid_xi, chebyshev_grid_yi, grid_xi[i], n);
		error[i] = std::abs(f(grid_xi[i]) - result_yi[i]);
	}
}

void murlib::build_least_sqr_sys(const int point_n, const int pol_deg, const double* xi, const double* yi, double* A, double* b) {
	// надо найти пседорешенеи Ах = b
	// множество псевдорешений совпадает с множеством решений А^T * A x = A^T b
	for (int i = 0; i < pol_deg; ++i) {
		for (int j = 0; j < pol_deg; ++j) {
			A[i * pol_deg + j] = 0;
			b[i] = 0;
			for (int k = 0; k < point_n; ++k) {
				A[i * pol_deg + j] += pow(xi[k], j + i);
				b[i] += pow(xi[k], i) * yi[k];
			}
		}
	}
}

double murlib::least_sqr_approx(const int deg, const double* coef, const double x) {
	double result = 0;
	for (int i = 0; i < deg; ++i) {
		result += pow(x, i) * coef[i];
	}
	return result;
}