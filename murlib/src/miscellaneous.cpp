
#include <cmath>
#include <iostream>
#include <functional>

#include "murlib/miscellaneous.h"
#include "murlib/matrix_solver.h"

double murlib::get_chebyshev_root(const double start, const double end, const int degree, const int i) {
	return (end + start) / 2 + (end - start) / 2 * cos(murlib::PI * (2 * i + 1) / 2 / degree);
}

void murlib::build_grid(double start, double end, double* test_xi, int grid_size) {
	for (int i = 0; i < grid_size; ++i) {
		test_xi[i] = start + (end - start) / (grid_size) * i;
	}
}

double murlib::determinant(const int n, double* matrix) {
	double* L = new double[n * n];
	double* U = new double[n * n];

	murlib::lu_decompostion(n, matrix, L, U);
	double result = 1;
	for (int i = 0; i < n; ++i) {
		result *= L[i + i * n] * U[i + i * n];
	}

	delete[] L;
	delete[] U;
	return result;
}

