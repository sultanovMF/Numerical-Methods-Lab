#pragma once

#include <iostream> 
#include <functional>
#include <cmath>
#include <tuple>
#include <iomanip>
#include <span>
#include <murlib/matrix_solver.h>
namespace murlib {
    std::pair<double, int> bisection(double a, double b, std::function<double(double)> f, const double eps = 1e-6) {
        if (f(a) * f(b) >= 0) {
            throw std::invalid_argument("You have not assumed right a and b");
        }
        int steps = 0;
        double c = a;
        while ((b - a) >= eps) {
            c = (a + b) / 2;
            if (f(c) * f(a) < 0)
                b = c;
            else
                a = c;

            ++steps;
        }

        return { c, steps };
    }

    std::pair<double, int> chord(double a, double b, std::function<double(double)> f, double eps = 1e-6) {
        int steps = 0;
        while (fabs(b - a) > eps) {
            a = a - (b - a) * f(a) / (f(b) - f(a));
            b = b - (a - b) * f(b) / (f(a) - f(b));
            steps++;
        }

        return { b, steps };
    }

    std::pair<double, int> simple_iter(double x_prev, std::function<double(double)> phi, double eps = 1e-6) {
        int steps = 0;
        double x = x_prev;
        do {
            x_prev = x;
            x = phi(x_prev);
            steps++;
        } while (abs(x - x_prev) >= eps);

        return { x, steps };
    }

    std::pair<double, int> secant(double x_prev, double x_prev_prev, std::function<double(double)> f, double eps = 1e-6) {
        int steps = 0;
        double x = x_prev;
        do {
            x_prev = x;
            x = x_prev - (x_prev - x_prev_prev) / (f(x_prev) - f(x_prev_prev)) * f(x_prev);
            x_prev_prev = x_prev;
            steps++;
        } while (abs(x - x_prev) >= eps);

        return { x, steps };
    }

    std::tuple<double, double, int> simple_iter(double x, double y, std::function<double(double, double)> psix, std::function<double(double, double)> psiy, double eps = 1e-6) {
        int steps = 0;
        double x_prev;
        double y_prev;
        do {
            x_prev = x;
            y_prev = y;

            x = psix(x_prev, y_prev);
            y = psiy(x_prev, y_prev);

            steps++;
        } while (std::max(abs(x - x_prev), abs(y - y_prev)) >= eps);

        return { x, y, steps };
    }


    std::tuple<double, double, int> newton_iter_exact(double x, double y,
        std::function<double(double, double)> f,
        std::function<double(double, double)> g,
        std::function<double(double, double)> fx,
        std::function<double(double, double)> fy,
        std::function<double(double, double)> gx,
        std::function<double(double, double)> gy,

        double eps = 1e-6) {
        int steps = 0;
        double x_prev;
        double y_prev;
        do {
            x_prev = x;
            y_prev = y;

            double matrix[4] = { fx(x_prev, y_prev), fy(x_prev, y_prev), gx(x_prev, y_prev), gy(x_prev, y_prev) };
            double b[2] = { -f(x_prev, y_prev), -g(x_prev, y_prev) };
            double delta[2];

            murlib::solve_gauss(2, matrix, b, delta);

            x = x_prev + delta[0];
            y = y_prev + delta[1];

            steps++;
        } while (std::max(abs(x - x_prev), abs(y - y_prev)) >= eps);

        return { x, y, steps };
    }

    std::tuple<double, double, int> newton_iter_approx(double x, double y,
        std::function<double(double, double)> f,
        std::function<double(double, double)> g,

        double eps = 1e-6) {
        int steps = 0;
        double x_prev;
        double y_prev;
        do {
            x_prev = x;
            y_prev = y;

            double matrix[4] = {
                (f(x_prev + eps, y_prev) - f(x_prev, y_prev)) / eps,
                (f(x_prev, y_prev + eps) - f(x_prev, y_prev)) / eps,
                (g(x_prev + eps, y_prev) - g(x_prev, y_prev)) / eps,
                (g(x_prev, y_prev + eps) - g(x_prev, y_prev)) / eps };

            double b[2] = { -f(x_prev, y_prev), -g(x_prev, y_prev) };
            double delta[2];

            murlib::solve_gauss(2, matrix, b, delta);

            x = x_prev + delta[0];
            y = y_prev + delta[1];

            steps++;
        } while (std::max(abs(x - x_prev), abs(y - y_prev)) >= eps);

        return { x, y, steps };
    }


}