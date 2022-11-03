
#include <iostream>
#include <functional>
#include <cmath>
// Библиотека для рандома
#include <effolkronium/random.hpp>

#include "murlib/interpolation.h"
#include "murlib/approximation.h"
#include "murlib/matrix_solver.h"
#include "murlib/integration.h"
#include "murlib/miscellaneous.h"

double murlib::integrate_rectangle(const double* x, const std::function<double(double)> f, const int n) {
	double result = 0;
	for (int i = 0; i < n - 1; ++i) {
		result +=f(x[i]) * (x[i + 1] - x[i]);
	}
	return result;
}

double murlib::integrate_trapezoid(const double* x, const std::function<double(double)> f, const int n) {
	double result = 0;
	for (int i = 0; i < n - 1; ++i) {
		result += (f(x[i]) + f(x[i + 1])) / 2 * (x[i + 1] - x[i]);
	}
	return result;
}

double murlib::integrate_simpson(const double* x, const std::function<double(double)> f, const int n)
{
	double h = x[1] - x[0];
	double result = 0;
	for (int i = 1; i < n; i+= 2) {
		result += (f(x[i - 1]) + 4 * f(x[i]) + f(x[i + 1]));
	}

	// result += (x[i + 1] - x[i]) / 6. * (f(x[i]) + 4. * f((x[i] + x[i + 1]) / 2.) + f(x[i + 1]));
	return result * h / 3.;
}



double murlib::aitken_trapezoid(double start, double end, const int n, const std::function<double(double)> f) {
	double* xh   = new double[n];
	double* xh2  = new double[2 * n];
	double* xh4 = new double[4 * n];

	murlib::build_grid(start, end, xh,  n);
	murlib::build_grid(start, end, xh2, 2 * n);
	murlib::build_grid(start, end, xh4, 4 * n);

	double B = (murlib::integrate_trapezoid(xh2, f, 2 * n) - murlib::integrate_trapezoid(xh, f, n))
		/ (murlib::integrate_trapezoid(xh4, f, 4 * n) - murlib::integrate_trapezoid(xh2, f, n * 2));

	delete[] xh;
	delete[] xh2;
	delete[] xh4;

	return std::log2(B);
}

double murlib::aitken_simpson(double start, double end, const int n, const std::function<double(double)> f) {
	double* xh = new double[n];
	double* xh2 = new double[2 * n];
	double* xh4 = new double[4 * n];

	murlib::build_grid(start, end, xh, n);
	murlib::build_grid(start, end, xh2, 2 * n);
	murlib::build_grid(start, end, xh4, 4 * n);

	double B = (murlib::integrate_simpson(xh2, f, 2 * n) - murlib::integrate_simpson(xh, f, n))
		/ (murlib::integrate_simpson(xh4, f, 4 * n) - murlib::integrate_simpson(xh2, f, n * 2));

	delete[] xh;
	delete[] xh2;
	delete[] xh4;

	return std::log2(B);
}


double murlib::aitken_rectangle(double start, double end, const int n, const std::function<double(double)> f) {
	double* xh = new double[n];
	double* xh2 = new double[2 * n];
	double* xh4 = new double[4 * n];

	murlib::build_grid(start, end, xh, n);
	murlib::build_grid(start, end, xh2, 2 * n);
	murlib::build_grid(start, end, xh4, 4 * n);

	double B = (murlib::integrate_rectangle(xh2, f, 2 * n) - murlib::integrate_rectangle(xh, f, n))
		/ (murlib::integrate_rectangle(xh4, f, 4 * n) - murlib::integrate_rectangle(xh2, f, n * 2));

	delete[] xh;
	delete[] xh2;
	delete[] xh4;

	return std::log2(B);
}

double murlib::runge_rectangle(double start, double end, const int n, const int p, const std::function<double(double)> f) {
	double sigma = pow(0.5, p) / (pow(0.5, p) - 1);
	double* xh = new double[n];
	double* xhk = new double[2 * n];

	murlib::build_grid(start, end, xh, n);
	murlib::build_grid(start, end, xhk, 2*n);
	double result = sigma * murlib::integrate_rectangle(xh, f, n) + (1 - sigma) * murlib::integrate_rectangle(xhk, f, 2 * n);

	delete[] xh;
	delete[] xhk;

	return result;
}

double murlib::romberg_rectangle(double start, double end, const int n, const int p, const int q, const std::function<double(double)> f) {
	double* matrix_first  = new double[q * q];
	double* matrix_second = new double[q * q];
	double* h = new double[q];

	for (int i = 0; i < q; ++i) {
		h[i] = (end - start) / (n * (i+1)) ;
	}

	for (int i = 0; i < q; ++i) {
		double* xi = new double[n * (i+1)];
		murlib::build_grid(start, end, xi, n * (i + 1));

		matrix_first[i * q] = murlib::integrate_rectangle(xi, f, n * (i + 1));
		matrix_second[i * q] = 1.;

		delete[] xi;
	}
	for (int i = 0; i < q; ++i) {
		for (int j = 1; j < q; ++j) {
			matrix_first [i * q + j] = pow(h[i], p + j - 1);
			matrix_second[i * q + j] = pow(h[i], p + j - 1);
		}
	}

	double result = determinant(q, matrix_first) / determinant(q, matrix_second);
	delete[] matrix_first;
	delete[] matrix_second;
	delete[] h;
	return result;
}

double murlib::adaptive_integration(double delta, double start, double end, const int n, double* xi, const std::function<double(double)> f) {
	double h = (end - start) / 2.;
	double S = 0.;
	
	xi[0] = start;
	
	for (int i = 1; i < n; ++i) {
		xi[i] = xi[i - 1] + h;
		while (true) {
			double Itr = h / 2 * (f(xi[i - 1]) + f(xi[i]));
			double Itr_sost = h / 4 * (f(xi[i - 1]) + 2 * f((xi[i - 1] + xi[i]) / 2.) + f(xi[i]));
			double eps = 1. / 3. * (Itr_sost - Itr);


			if (abs(eps) > h * delta / (end - start)) {
				h /= 2.;
				xi[i] = xi[i - 1] + h;
			}
			else {
				xi[i] = xi[i - 1] + h;
				S += Itr_sost;
				break;
			}
		}
		if (xi[i] + h > end) {
			h = end - xi[i];
		}
	}

	return S;
}

double murlib::monte_carlo(double start, double end, const int n, const std::function<double(double)> f) {
	double result = 0.;
	{
		using Random = effolkronium::random_static;
		for (int i = 0; i < n; ++i) {
			double randx = Random::get(start, end);
			result += f(randx);
		}
	}
	return (end - start) * result / n;

}