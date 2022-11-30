#pragma once

namespace murlib {
	void tridiagonal_matrix(const int n, const double* A, const double* B, const double* C, const double* D,  double* X);
	void pentadiagonal_matrix(
		const int n, 
		const double* A,
		const double* B,
		const double* C,
		const double* D,
		const double* E,
		const double* F,
		double* X);
	void solve_gauss(const unsigned int n, double* A, double* b, double* x);

	void solve_lu(const unsigned int n, double* L, double* U, double* b, double* x);
	void lu_decompostion(const unsigned int n, double* A, double* L, double* U);
	void sqroot_decompostion(const unsigned int n, double* A, double* S);
	void sqroot_solve(const unsigned int n, double* S, double* b, double* x);

	void transpose(const unsigned int n, double* A, double* AT);
	void inverse(const int n, double* A, double* result);
	void multiply(const int n, double* A, double* B, double* C);
	void multiply_to_vec(const int n, const double* A, const double* x, double* b);
}