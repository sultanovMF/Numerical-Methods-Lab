#include <functional>



namespace murlib {
	void mnrp_for_polynoms(const std::function<double(double)> f, const int n, const double a, const double b, const int test_grid_size, const double* grid_xi, double* result_yi, double* error);
	void build_least_sqr_sys(const int point_n, const int pol_deg, const double* xi, const double* yi, double* A, double* b);
	double least_sqr_approx(const int deg, const double* coef, const double x);
}

