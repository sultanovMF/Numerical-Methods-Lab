
#include <cmath>
#include <iostream>
#include <functional>

#include "murlib/miscellaneous.h"


double murlib::get_chebyshev_root(const double start, const double end, const int degree, const int i) {
	return (end + start) / 2 + (end - start) / 2 * cos(murlib::PI * (2 * i + 1) / 2 / degree);
}
