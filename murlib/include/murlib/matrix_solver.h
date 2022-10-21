#pragma once


namespace murlib {
	void tridiagonal_matrix(const int n, const double* A, const double* B, const double* C, const double* D,  double* X);
	void solve_gauss(const unsigned int n, double* A, double* b, double* x);

	void solve_lu(const unsigned int n, double* L, double* U, double* b, double* x);
	void lu_decompostion(const unsigned int n, double* A, double* L, double* U);


}