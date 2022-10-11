#include <cmath>
#include <functional>

namespace murlib {
	double integrate_rectangle(const double* x, const double* y, const int n);
	double integrate_trapezoid(const double* x, const double* y, const int n);
	double integrate_simpson(const double* x, const std::function<double(double)> f, const int n);
	double aitken_process(const double a, const double b, const double k, const double h, const std::function<double(double)> f);
	double aitken_process(const double a, const double b, const double k, const double h, const std::function<double(double)> f, double& S_1, double& S1);
	double runge_method(const double S_1, const double S_2, const double k, const double p);
	// TODO: Romberg
	//

}