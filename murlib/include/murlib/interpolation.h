#pragma once


namespace murlib {
    const double PI = 3.14159265358979323846;
    double lagrange_polynom(const double* xi, const double* yi, const double t, const int n);
    double trigonometric_polynom(double t, const double xi[], const double yi[], const int n);
    double newton_polynom(const double* xi, double** v, const double t, const int degree);
    double get_chebyshev_root(const double start, const double end, const int degree, const int i);

    void build_splines(const int n, double* A, double* B, double* C, double* D, const double* X, const double* Y);
    double spline_iterpolation(const int n, const double* A, const  double* B, const  double* C, const  double* D, const double* X, const int i, const double t);
}