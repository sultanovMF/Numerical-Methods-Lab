#pragma once

#include <vector>

namespace murlib {
	double lagrange_polynom(const std::vector<double>& x, const std::vector<double>& y, const double t, const int n);
    double newton_polynom(const std::vector<double>& x, const std::vector<double>& y, const double t, const int n);
}

double murlib::lagrange_polynom(const std::vector<double>& x, const std::vector<double>& y, const double t, const int n) {
    double result = 0;

    for (int i = 0; i < n; ++i)
    {
        double term = y[i];
        for (int j = 0; j < n; ++j)
        {
            if (j != i) {
                term = term * (t - x[j]) / (x[i] - x[j]);
            }
        }

        result += term;
    }

    return result;
}

double murlib::newton_polynom(const std::vector<double>& x, const std::vector<double>& y, const double t, const int n) {
    double result = y[0], F, den;
    for (int i = 1; i < n; ++i) {
        F = 0;
        for (int j = 0; j <= i; ++j) {
            den = 1;
            for (int k = 0; k <= i; k++) {
                if (k != j) den *= (x[j] - x[k]);
            }
            //считаем разделенную разность
            F += y[j] / den;
        }
        //домножаем разделенную разность на скобки (x-x[0])...(x-x[i-1])
        for (int k = 0; k < i; k++)F *= (t - x[k]);
        result += F;//полином
    }
    return result;
}
