#pragma once

#include <iostream>
#include <functional>

#include "murlib/interpolation.h"
#include "murlib/approximation.h"
#include "murlib/matrix_solver.h"
#include "murlib/integration.h"


double murlib::integrate_rectangle(const double* x, const double* y, const int n) {
	double result = 0;
	for (int i = 0; i < n - 1; ++i) {
		result += y[i] * (x[i + 1] - x[i]);
	}
	return result;
}

double murlib::integrate_trapezoid(const double* x, const double* y, const int n) {
	double result = 0;
	for (int i = 0; i < n - 1; ++i) {
		result += (y[i] + y[i + 1]) / 2 * (x[i + 1] - x[i]);
	}
	return result;
}

double murlib::integrate_simpson(const double* x, const std::function<double(double)> f, const int n)
{
	double result = 0;
	for (int i = 0; i < n - 1; ++i) {
		result += (f(x[i]) + 4. * f((x[i] + x[i + 1]) / 2.) + f(x[i + 1]) / 2. * (x[i + 1] - x[i]));
	}
	result /= 6.;
	return 0.0;
}



double murlib::aitken_process(const double a, const double b, const double k, const double h, const std::function<double(double)> f, double& S_1, double& S_2) {
	const int n_1 = std::ceil((a + b) / h);
	const int n_2 = std::ceil((a + b) / h / k);
	const int n_3 = std::ceil((a + b) / h / std::pow(k, 2));

	double* grid_xi_1 = new double[n_1];
	double* grid_xi_2 = new double[n_2];
	double* grid_xi_3 = new double[n_3];

	double* grid_yi_1 = new double[n_1];
	double* grid_yi_2 = new double[n_2];
	double* grid_yi_3 = new double[n_3];

	for (int i = 0; i < n_1; ++i) {
		grid_xi_1[i] = a + (b - a) / n_1 * i;
		grid_yi_1[i] = f(grid_xi_1[i]);
	}
	for (int i = 0; i < n_2; ++i) {
		grid_xi_2[i] = a + (b - a) / n_2 * i;
		grid_yi_2[i] = f(grid_xi_2[i]);
	}
	for (int i = 0; i < n_3; ++i) {
		grid_xi_3[i] = a + (b - a) / n_3 * i;
		grid_yi_3[i] = f(grid_xi_3[i]);
	}

	const double S = murlib::integrate_simpson(grid_xi_1, f, n_1);
	S_1 = murlib::integrate_simpson(grid_xi_2, f, n_2);
	S_2 = murlib::integrate_simpson(grid_xi_3, f, n_3);

	const double B = (S_1 - S_2) / (S_2 - S_2);

	delete[] grid_xi_1;
	delete[] grid_xi_2;
	delete[] grid_xi_3;

	delete[] grid_yi_1;
	delete[] grid_yi_2;
	delete[] grid_yi_3;

	return -log(B) / log(k);
}
double murlib::aitken_process(const double a, const double b, const double k, const double h, const std::function<double(double)> f) {
	double S_1, S_2;
	return murlib::aitken_process(a, b, k, h, f, S_1, S_2);
}


double murlib::runge_method(const double S_1, const double S_2, const double k, const double p) {
	const double sigma = pow(k, p) / (pow(k, p) - 1);
	return sigma * S_1 + (1 - sigma) * S_2;
}

