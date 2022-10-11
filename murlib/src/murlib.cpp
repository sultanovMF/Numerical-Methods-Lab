#pragma once

#include <cmath>
#include <iostream>
#include <functional>

#include "murlib/interpolation.h"
#include "murlib/approximation.h"
#include "murlib/matrix_solver.h"
#include "murlib/integration.h"



double murlib::lagrange_polynom(const double* xi, const double* yi, const double t, const int n) {
	double result = 0;
	for (int i = 0; i < n; ++i)
	{
		double term = yi[i];
		for (int j = 0; j < n; ++j)
		{
			if (j != i) {
				term = term * (t - xi[j]) / (xi[i] - xi[j]);
			}
		}
		result += term;
	}
	return result;
}

double murlib::newton_polynom(const double* xi, double** v, const double t, const int n) {
	// А я же каждый раз пересчитываю v получается...
	// TODO Надо расчет v вынести в отдельную функцию
	// TODO Добавить возможность добавлять узлы интерполяции
	for (int i = 1; i < n; ++i) {
		for (int j = 0; j < n - i; ++j) {
			v[i][j] = (v[i - 1][j + 1] - v[i - 1][j]) / (xi[j + i] - xi[j]);
		}
	}

	double result = v[0][0];
	double term = 1;
	for (int i = 1; i < n; ++i) {
		term *= (t - xi[i - 1]);
		result += v[i][0] * term;
	}

	return result;
}

double murlib::trigonometric_polynom(double x, const double xi[], const double yi[], const int n) {
	// n - степень многочлена
	double points_number = 2 * n + 1;
	
	// Инциализация начальных данных
	double a_0 = 0.;
	double* a = new double[n];
	double* b = new double[n];
	for (int i = 0; i < n; ++i) {
		a[i] = 0;
		b[i] = 0;
	}
	//СЛЕДУЮЩИЙ КОД РАБОЧИЙ!!
	for (int m = 0; m < n; ++m) {
		for (int k = 0; k < points_number; ++k) {
			if (m == 0) a_0 += yi[k];
			a[m] += yi[k] * cos((m + 1) * xi[k]);
			b[m] += yi[k] * sin((m + 1) * xi[k]);
		}
		if (m == 0) a_0 = a_0 * 2 / (points_number);
		a[m] = a[m] * 2 / (points_number);
		b[m] = b[m] * 2 / (points_number);
	}
	double result = a_0 / 2.;
	for (int k = 0; k < n; ++k) {
		result += a[k] * cos(x * (k+1)) + b[k] * sin(x * (k+1));
	}

	delete[] a;
	delete[] b;

	return result;
}

double murlib::get_chebyshev_root(const double start, const double end, const int degree, const int i) {
	return (end + start) / 2 + (end - start) / 2 * cos(murlib::PI * (2 * i + 1) / 2 / degree);
}


void murlib::tridiagonal_matrix(const int n, const double* A, const double* B, const double* C, const double* D, double* X) {
	// Еще есть алгоритм быстрее где можно изменять А б с
	double* P = new double[n-1];
	double* Q = new double[n];

	// forward
	P[0] = C[0] / B[0];
	Q[0] = D[0] / B[0];
	for (int i = 1; i < n; ++i) {
		if (i < n-1) P[i] = C[i] / (B[i] - A[i] * P[i - 1]);
		Q[i] = (D[i] - A[i] * Q[i-1]) / (B[i] - A[i] * P[i - 1]);
	}

	//backward
	X[n - 1] = Q[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		X[i] = Q[i] - P[i] * X[i + 1];
	}
	delete[] P;
	delete[] Q;
}

double murlib::spline_iterpolation(const int n, const double* A, const  double* B, const  double* C, const  double* D, const double* X, const int i, const double t) {
	return A[i] + B[i] * (t - X[i]) + C[i] * (t - X[i]) * (t - X[i]) + D[i] * (t - X[i]) * (t - X[i]) * (t - X[i]);
}

void murlib::build_splines(const int n, double* A, double* B, double* C, double* D, const double X[], const double Y[]) {
	// C должен быть размера n

	double* H = new double[n];
	for (int i = 0; i < n; ++i) {
		H[i] = X[i + 1] - X[i];
	}
	// если n = 2, то уравнений 8, то есть система будет размера 4N неизвестны

	double* AA = new double[n + 1];
	double* BB = new double[n + 1];
	double* CC = new double[n + 1];
	double* DD = new double[n + 1];

	AA[0] = 0;
	BB[0] = 1;
	CC[0] = 0;
	DD[0] = 0;
	for (int i = 1; i < n; ++i) {
		BB[i] = 2 * (H[i - 1] + H[i]);
		CC[i] = H[i];
		DD[i] = 3 * (Y[i + 1] - Y[i]) / H[i] - 3 * (Y[i] - Y[i - 1]) / H[i - 1];
		AA[i] = H[i-1];
	}
	AA[n] = 0;
	CC[n] = 0;
	DD[n] = 0;

	double* CCC = new double[n + 1]; // то куда мы положим результат трехдиагональной прогонки
	murlib::tridiagonal_matrix(n + 1, AA, BB, CC, DD, CCC);
	
	for (int i = 0; i < n; ++i) {
		A[i] = Y[i];
		B[i] = ((Y[i+1] - Y[i]) / H[i]) - H[i] * (CCC[i + 1] + 2 * CCC[i]) / 3.;
		C[i] = CCC[i];
		D[i] = (CCC[i + 1] - CCC[i]) / 3. / H[i];
	}

	delete[] H;
	delete[] AA;
	delete[] BB;
	delete[] CC;
	delete[] DD;
	delete[] CCC;
}

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
		result += (f(x[i]) + 4. * f((x[i] + x[i+1]) /2.) + f(x[i+1]) / 2. * (x[i + 1] - x[i]));
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

