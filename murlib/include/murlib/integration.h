#include <cmath>
#include <functional>

namespace murlib {
	double integrate_rectangle(const double* x, const std::function<double(double)> f, const int n);
	double integrate_trapezoid(const double* x, const std::function<double(double)> f, const int n);
	double integrate_simpson(const double* x, const std::function<double(double)> f, const int n);

	double aitken_trapezoid (double start, double end, const int n, const std::function<double(double)> f);
	double aitken_rectangle (double start, double end, const int n, const std::function<double(double)> f);
	double aitken_simpson   (double start, double end, const int n, const std::function<double(double)> f);
	double runge_rectangle  (double start, double end, const int n, const int p, const std::function<double(double)> f);
	double romberg_rectangle(double start, double end, const int n, const int p, const int q, const std::function<double(double)> f);
	double adaptive_integration(double delta, double start, double end, const int n, double* xi, const std::function<double(double)> f);
	double gauss_quad(double N, double m, const std::function<double(double)> f);
	double monte_carlo(double start, double end, const int n, const std::function<double(double)> f);
}