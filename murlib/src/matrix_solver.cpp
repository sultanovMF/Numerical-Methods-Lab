
#include <cmath>
#include <stdexcept>

#include "murlib/matrix_solver.h"
#include "murlib/miscellaneous.h"

void backward_up(const unsigned int n, double* U, double* b, double* x) {
	for (int i = n - 1; i >= 0; --i) {
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= U[j + i * n] * x[j];
		}
		x[i] /= U[i + i * n];
	}
}

void backward_low(const unsigned int n, double* L, double* b, double* y) {
	for (int i = 0; i < n; ++i) {
		y[i] = b[i];
		for (int j = 0; j < i; ++j) {
			y[i] -= L[i * n + j] * y[j];
		}
		y[i] /= L[i + i * n];
	}
}

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
		x[i] = b[i];
		for (int j = i + 1; j < n; ++j) {
			x[i] -= A[j + i * n] * x[j];
		}
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

			L[j * n + i] = A[j * n + i];

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
	backward_low(n, L, b, y);
	backward_up(n, U, y, x);
	delete[] y;
}
void murlib::transpose(const unsigned int n, double* A, double* AT) {
	std::copy(A, A + n*n, AT);
	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j < n; ++j) {
			std::swap(AT[i * n + j], AT[j * n + i]);
		}
	}
}

void murlib::sqroot_decompostion(const unsigned int n, double* A, double* S) {
	std::fill(S, S + n * n, 0.);
	for (int i = 0; i < n; ++i) {
		S[i + i * n] = A[i + i * n];

		for (int p = 0; p < i; ++p) {
			S[i + i * n] -= S[i * n + p] * S[i * n + p];
		}

		S[i + i * n] = std::sqrt(S[i + i * n]);

		for (int j = i + 1; j < n; ++j) {
			S[j * n + i] = A[j * n + i];

			for (int p = 0; p < i; ++p) {
				S[j * n + i] -= S[i * n + p] * S[j * n + p];
			}

			S[j * n + i] /= S[i + i * n];
		}
	}
}

void murlib::sqroot_solve(const unsigned int n, double* S, double*b, double* x) {
	double* ST = new double[n * n];
	transpose(n, S, ST);

	double* y = new double[n];
	backward_low(n, S, b, y); //-59 -33 -12


	backward_up(n, ST, y, x);
	delete[] y;
	delete[] ST;
}


void murlib::pentadiagonal_matrix(
	const int n,
	const double* A,
	const double* B,
	const double* C,
	const double* D,
	const double* E,
	const double* F,
	double* x) {

	double* alpha = new double	[n];
	double* beta = new double	[n];
	double* gamma = new double	[n];
	double z;

	z = C[0];
	alpha[0] = -D[0] / z;
	beta[0]	 = -E[0] / z;
	gamma[0] =  F[0] / z;

	z = B[1] * alpha[0] + C[1];
	alpha[1] = -(B[1] * beta[0] + D[1]) / z;
	beta[1]  = -(E[1] / z);
	gamma[1] = (F[1] - B[1] * gamma[0]) / z;

	for (int i = 2; i < n; ++i) {
		// Вообще говоря, здесь идет лишний счет, TODO переделай!!!
		z = A[i] * alpha[i - 2] * alpha[i - 1] + A[i] * beta[i - 2] + B[i] * alpha[i - 1] + C[i];
		alpha[i] = -(A[i] * alpha[i - 2] * beta[i - 1] + B[i] * beta[i - 1] + D[i]) / z;
		beta[i]  = -(E[i]) / z;
		gamma[i] = (F[i] - (A[i] * alpha[i - 2] * gamma[i - 1] + A[i] * gamma[i - 2] + B[i] * gamma[i - 1])) / z;
	}

	z = (A[n - 1] * alpha[n - 3] * alpha[n - 2] + A[n - 1] * beta[n - 3] + B[n - 1] * alpha[n - 2] + C[n - 1]);
	
	x[n - 1] = (F[n - 1] - (A[n - 1] * alpha[n - 3] * gamma[n - 2] + A[n-1] * gamma[n - 3] + B[n - 1] * gamma[n - 2]))
		/ (A[n - 1] * alpha[n - 3] * alpha[n - 2] + A[n - 1] * beta[n - 3] + B[n - 1] * alpha[n - 2] + C[n - 1]);
	// n = 6
	x[n - 2] = alpha[n - 2] * x[n - 1] + gamma[n - 2];

	for (int i = n - 3; i >= 0; --i) {
		x[i] = alpha[i] * x[i + 1] + beta[i] * x[i + 2] + gamma[i];
	}

	delete[] alpha;
	delete[] beta;
	delete[] gamma;
}

void murlib::inverse(const int n, double* A, double* result) {
	// TODO testing!
	double* L = new double[n * n];
	double* U = new double[n * n];

	murlib::lu_decompostion(n, A, L, U);

	for (int i = 0; i < n; ++i) {
		double* b = new double[n];
		std::fill(b, b + n, 0);
		b[i] = 1;

		double* x = new double[n];
		
		murlib::solve_lu(n, L, U, b, x);

		for (int j = 0; j < n; ++j) {
			result[j * n + i] = x[j];
		}

		delete[] b;
	}

	delete[] L;
	delete[] U;
}

void murlib::multiply_to_vec(const int n, const double* A, const double* x, double* b) {
	for (int i = 0; i < n; ++i) {
		b[i] = 0.;
		for (int j = 0; j < n; ++j) {
			b[i] += A[j + i * n] * x[j];
		}
	}
}

void murlib::multiply(const int N, double* A, double* B, double* C) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < N; ++k) {
				C[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}
}


///// N - размерность матрицы; A[N][N] - матрица коэффициентов, F[N] - столбец свободных членов,
///// X[N] - начальное приближение, ответ записывается также в X[N];
//void Jacobi(int N, double** A, double* F, double* X)
//{
//	double* TempX = new double[N];
//	double norm; // норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.
//
//	do {
//		for (int i = 0; i < N; i++) {
//			TempX[i] = F[i];
//			for (int g = 0; g < N; g++) {
//				if (i != g)
//					TempX[i] -= A[i][g] * X[g];
//			}
//			TempX[i] /= A[i][i];
//		}
//		norm = fabs(X[0] - TempX[0]);
//		for (int h = 0; h < N; h++) {
//			if (fabs(X[h] - TempX[h]) > norm)
//				norm = fabs(X[h] - TempX[h]);
//			X[h] = TempX[h];
//		}
//	} while (norm > eps);
//	delete[] TempX;
//}


//Inputs: A, b, ω
//Output : φ
//
//Choose an initial guess φ to the solution
//repeat until convergence
//	for i from 1 until n do
//		set σ to 0
//		for j from 1 until n do
//			if j ≠ i then
//				set σ to σ + aij φj
//			end if
//		end(j - loop)
//		set φi to(1 − ω)φi + ω(bi − σ) / aii
//	end(i - loop)
//	check if convergence is reached
//end(repeat)