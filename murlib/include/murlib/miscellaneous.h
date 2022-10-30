#pragma once 

namespace murlib {
	constexpr double PI = 3.14159265358979323846;
	constexpr double EPSILON = 0.0000001;

	double get_chebyshev_root(const double start, const double end, const int degree, const int i);
}