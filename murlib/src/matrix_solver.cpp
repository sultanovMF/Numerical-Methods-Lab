#pragma once

#include <cmath>
#include <stdexcept>

#include "murlib/matrix_solver.h"
#include "murlib/miscellaneous.h"

void murlib::tridiagonal_matrix(const int n, const double* A, const double* B, const double* C, const double* D, double* X) {
	// Еще есть алгоритм быстрее где можно изменять А б с
	double* P = new double[n - 1];
	double* Q = new double[n];

	// forward
	P[0] = C[0] / B[0];
	Q[0] = D[0] / B[0];
	for (int i = 1; i < n; ++i) {
		if (i < n - 1) P[i] = C[i] / (B[i] - A[i] * P[i - 1]);
		Q[i] = (D[i] - A[i] * Q[i - 1]) / (B[i] - A[i] * P[i - 1]);
	}

	//backward
	X[n - 1] = Q[n - 1];
	for (int i = n - 2; i >= 0; --i) {
		X[i] = Q[i] - P[i] * X[i + 1];
	}
	delete[] P;
	delete[] Q;
}

void murlib::solve_gauss(const unsigned int n, double* A, double* b, double* x) {
	// не учитывает что могут быть пустые строки в матрице
	// forward
	for (int i = 0; i < n; ++i) {
		double pivot = A[i + i * n];
		if (abs(pivot) < murlib::EPSILON) {
			// Если пилотный элемент равен нулю
			double max = 0.;
			int max_index = i;
			for (int j = i + 1; j < n; ++j) {
				if (abs(A[i + j * n]) > max) {
					max_index = j;
					max = A[j + j * n];
				}
			}

			if (max < murlib::EPSILON) {
				throw std::invalid_argument("det(A) = 0");
			}
			for (int k = i; k < n; ++k) {
				std::swap(A[k + max_index * n], A[k + i * n]);
			}
			std::swap(b[i], b[max_index]);
			i--;
			continue;
		}
		for (int j = i; j < n; ++j) {
			A[j + i * n] /= pivot;
		}
		b[i] /= pivot;

		for (int j = i + 1; j < n; ++j) {
			pivot = A[i + j * n];
			for (int k = 0; k < n; ++k) {
				A[k + j * n] -= pivot * A[k + n * i];
			}
			b[j] -= pivot * b[i];
		}
	}

	//backward
	for (int i = n - 1; i >= 0; --i) {
		x[i] = 0;
		for (int j = i + 1; j < n; ++j) {
			x[i] -= A[j + i * n] * x[j];
		}
		x[i] += b[i];
	}
}

void murlib::lu_decompostion(const unsigned int n, double* A, double* L, double* U) {
	for (int i = 0; i < n * n; ++i) {
		L[i] = 0;
		U[i] = 0;
	}

	for (int i = 0; i < n; ++i) {
		L[i * n] = A[i * n];
		U[i] = A[i] / L[0];
	}

	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = A[i * n + j];

			for (int k = 0; k < i; k++)
			{
				U[i * n + j] -= L[i *n + k] * U[k * n + j];
			}

			L[j * n + i] = A[j * n + i];/* getSumL(n, i, j, mL, mU)) / mU[i][i];*/

			for (int k = 0; k < i; k++)
			{
				L[j * n + i] -= L[j*n+k] * U[k * n + i];
			}

			L[j * n + i] /= U[i * n + i];
		}
	}

}

void murlib::solve_lu(const unsigned int n, double* L, double* U, double* b, double* x) {
	double* y = new double[n];

	for (int i = 0; i < n; ++i) {
		y[i] = b[i];
		for (int j = 0; j < i; ++j) {
			y[i] -= L[i * n + j] * y[j];
		}
		y[i] /= L[i + i * n];
	}

	//for (int i = n - 1; i >= 0; --i) {
	//	x[i] = y[i];

	//	for (int j = i + 1; j < n; ++j) {
	//		x[i] -= U[i * n + j] * x[j];
	//	}


	//}
	for (int i = 0; i < n; ++i) {
		double pivot = U[i * n + i];
		for (int j = i; j < n; ++j) {
			U[i * n + j] /= pivot;
		}
		y[i] /= pivot;
	}

	for (int i = n - 1; i >= 0; --i) {
		x[i] = y[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= U[j + i * n] * x[j];
		}
	}
	//1 2 3
	//0 -2 -5
	//0 0 -1
	delete[] y;
}