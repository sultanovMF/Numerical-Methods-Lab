#pragma once

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
		error[i] = abs(f(grid_xi[i]) - result_yi[i]);
	}
}