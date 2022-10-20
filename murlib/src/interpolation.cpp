#pragma once

#include <cmath>
#include <iostream>
#include <functional>

#include "murlib/interpolation.h"
#include "murlib/approximation.h"
#include "murlib/matrix_solver.h"
#include "murlib/integration.h"
#include "murlib/miscellaneous.h"

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
		result += a[k] * cos(x * (k + 1)) + b[k] * sin(x * (k + 1));
	}

	delete[] a;
	delete[] b;

	return result;
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
		AA[i] = H[i - 1];
	}
	AA[n] = 0;
	CC[n] = 0;
	DD[n] = 0;

	double* CCC = new double[n + 1]; // то куда мы положим результат трехдиагональной прогонки
	murlib::tridiagonal_matrix(n + 1, AA, BB, CC, DD, CCC);

	for (int i = 0; i < n; ++i) {
		A[i] = Y[i];
		B[i] = ((Y[i + 1] - Y[i]) / H[i]) - H[i] * (CCC[i + 1] + 2 * CCC[i]) / 3.;
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