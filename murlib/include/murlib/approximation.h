#include <functional>



namespace murlib {
	void mnrp_for_polynoms(const std::function<double(double)> f, const int n, const double a, const double b, const int test_grid_size, const double* grid_xi, double* result_yi, double* error);
}

