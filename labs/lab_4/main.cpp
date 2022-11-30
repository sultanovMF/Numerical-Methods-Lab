#include <iostream> 
#include <functional>
#include <cmath>
#include <tuple>
#include <iomanip>
#include <span>
#include <murlib/matrix_solver.h>
#include <murlib/nonlinear_eq.h>

int main() {
    std::cout << std::setprecision(7);

    auto f = [](double x) {
        return 2. * std::sin(x * x) - 1.;
    };

    // Метод бисекций 
    {
        double a;
        double b;
        double c;
        int steps;
        // #1
        a = 4.4;
        b = 4.5;

        std::tie(c, steps) = murlib::bisection(a, b, f);
        std::cout << "Bisection root #1: " << c << " steps: " << steps << std::endl;

        // #2
        a = 4.6;
        b = 4.7;

        std::tie(c, steps) = murlib::bisection(a, b, f);
        std::cout << "Bisection root #2: " << c << " steps: " << steps << std::endl;

        // #3
        a = 5;
        b = 5.1;

        std::tie(c, steps) = murlib::bisection(a, b, f);
        std::cout << "Bisection root #3: " << c << " steps: " << steps << std::endl;
        // #4
        a = 5.2;
        b = 5.3;

        std::tie(c, steps) = murlib::bisection(a, b, f);
        std::cout << "Bisection root #4: " << c << " steps: " << steps << std::endl;

        // #5
        a = 5.6;
        b = 5.7;

        std::tie(c, steps) = murlib::bisection(a, b, f);
        std::cout << "Bisection root #5: " << c << " steps: " << steps << std::endl;

        // #6
        a = 5.8;
        b = 5.9;

        std::tie(c, steps) = murlib::bisection(a, b, f);
        std::cout << "Bisection root #6: " << c << " steps: " << steps << std::endl;
    }
    std::cout << std::endl;
    // Метод хорд
    {

        double a;
        double b;
        double c;
        int steps;
        // #1
        a = 4.4;
        b = 4.5;

        std::tie(c, steps) = murlib::chord(a, b, f);
        std::cout << "Chord root #1: " << c << " steps: " << steps << std::endl;

        // #2
        a = 4.6;
        b = 4.7;

        std::tie(c, steps) = murlib::chord(a, b, f);
        std::cout << "Chord root #2: " << c << " steps: " << steps << std::endl;

        // #3
        a = 5;
        b = 5.1;

        std::tie(c, steps) = murlib::chord(a, b, f);
        std::cout << "Chord root #3: " << c << " steps: " << steps << std::endl;
        // #4
        a = 5.2;
        b = 5.3;

        std::tie(c, steps) = murlib::chord(a, b, f);
        std::cout << "Chord root #4: " << c << " steps : " << steps << std::endl;

        // #5
        a = 5.6;
        b = 5.7;

        std::tie(c, steps) = murlib::chord(a, b, f);
        std::cout << "Chord root #5: " << c << " steps: " << steps << std::endl;

        // #6
        a = 5.8;
        b = 5.9;

        std::tie(c, steps) = murlib::chord(a, b, f);
        std::cout << "Chord root #6: " << c << " steps: " << steps << std::endl;
    }
    std::cout << std::endl;
    // Метод простых итераций
    {
        double x;
        int steps;
        {
            // Положим alpha(x) = const = 1 / 100
            // Тогда phi = 0.02 sin(x^2) - 0.01 + x
            // |phi'| < 1 при x = 4.6..4.7
            //                x = 5.2..5.3
            //                x = 5.8..5.9
            auto phi = [](double x) {
                return 0.02 * std::sin(x * x) - 0.01 + x;
            };
            x = 4.65;
            std::tie(x, steps) = murlib::simple_iter(x, phi);
            std::cout << "Simple iteration #2 " << x << " steps: " << steps << std::endl;

            x = 5.2;
            std::tie(x, steps) = murlib::simple_iter(x, phi);
            std::cout << "Simple iteration #4 " << x << " steps: " << steps << std::endl;

            x = 5.8;
            std::tie(x, steps) = murlib::simple_iter(x, phi);
            std::cout << "Simple iteration #6 " << x << " steps: " << steps << std::endl;
        }
        {
            // Положим alpha(x) = const = -1 / 100
            // Тогда phi = 0.02 sin(x^2) - 0.01 + x
            // |phi'| < 1 при x = 4.3..4.5
            //                x = 5.0..5.1
            //                x = 5.6..5.7

            auto phi = [](double x) {
                return x - 0.02 * std::sin(x * x) + 0.01;
            };

            x = 4.3;
            std::tie(x, steps) = murlib::simple_iter(x, phi);
            std::cout << "Simple iteration #1 " << x << " steps: " << steps << std::endl;

            x = 5.0;
            std::tie(x, steps) = murlib::simple_iter(x, phi);
            std::cout << "Simple iteration #3 " << x << " steps: " << steps << std::endl;

            x = 5.7;
            std::tie(x, steps) = murlib::simple_iter(x, phi);
            std::cout << "Simple iteration #5 " << x << " steps: " << steps << std::endl;
        }
    }
   
    std::cout << std::endl;
    // Метод Ньютона
    {
        double x;
        int steps;

        auto phi = [](double x) {
            // phi = x + f'(x)
            return x - (2 * sin(x * x) - 1) / (4. * x * std::cos(x * x));
        };

        x = 4.4;
        std::tie(x, steps) = murlib::simple_iter(x, phi);
        std::cout << "Newton iteration #1 " << x << " steps: " << steps << std::endl;

        x = 4.65;
        std::tie(x, steps) = murlib::simple_iter(x, phi);
        std::cout << "Newton iteration #2 " << x << " steps: " << steps << std::endl;

        x = 5.0;
        std::tie(x, steps) = murlib::simple_iter(x, phi);
        std::cout << "Newton iteration #3 " << x << " steps: " << steps << std::endl;

        x = 5.2;
        std::tie(x, steps) = murlib::simple_iter(x, phi);
        std::cout << "Newton iteration #4 " << x << " steps: " << steps << std::endl;

        x = 5.7;
        std::tie(x, steps) = murlib::simple_iter(x, phi);
        std::cout << "Newton iteration #5 " << x << " steps: " << steps << std::endl;

        x = 5.8;
        std::tie(x, steps) = murlib::simple_iter(x, phi);
        std::cout << "Newton iteration #6 " << x << " steps: " << steps << std::endl;
    }
    std::cout << std::endl;
    // Метод секущих
    {
        double x;
        double x_prev, x_prev_prev;
        int steps;


        x_prev = 4.3;
        x_prev_prev = 4.4;
        std::tie(x, steps) = murlib::secant(x_prev, x_prev_prev, f);
        std::cout << "Secant iteration #1 " << x << " steps: " << steps << std::endl;

        x_prev = 4.64;
        x_prev_prev = 4.65;
        std::tie(x, steps) = murlib::secant(x_prev, x_prev_prev, f);
        std::cout << "Secant iteration #2 " << x << " steps: " << steps << std::endl;

        x_prev = 5.0;
        x_prev_prev = 5.05;
        std::tie(x, steps) = murlib::secant(x_prev, x_prev_prev, f);
        std::cout << "Secant iteration #3 " << x << " steps: " << steps << std::endl;

        x_prev = 5.2;
        x_prev_prev = 5.23;
        std::tie(x, steps) = murlib::secant(x_prev, x_prev_prev, f);
        std::cout << "Secant iteration #4 " << x << " steps: " << steps << std::endl;

        x_prev = 5.7;
        x_prev_prev = 5.75;
        std::tie(x, steps) = murlib::secant(x_prev, x_prev_prev, f);
        std::cout << "Secant iteration #5 " << x << " steps: " << steps << std::endl;

        x_prev = 5.78;
        x_prev_prev = 5.8;
        std::tie(x, steps) = murlib::secant(x_prev, x_prev_prev, f);
        std::cout << "Secant iteration #6 " << x << " steps: " << steps << std::endl;
    }
    // МПИ для поиска минимума
    {
        // Условия минимума:
        // psi_x = 0
        // psi_y = 0
        auto psix = [](double x, double y) {
            return  x - (2 * x + tan(x + y) + (x + y + 1) * (1 + tan(x + y) * tan(x + y))) / 10'000;
        };

        auto psiy = [](double x, double y) {
            return y - (2 * y + tan(x + y) + (x + y + 1) * (1 + tan(x + y) * tan(x + y))) / 10'000;
        };
        double x = -0.176;
        double y = -0.176;
        int steps;
        std::tie(x, y, steps) = murlib::simple_iter(x, y, psix, psiy);
        std::cout << "Min simple iterations x: " << x << " y: " << y << " steps: " << steps << std::endl;
    }
    std::cout << std::endl;
    // Метод Ньютона для поиска минимума
    {  
        double x = -0.176;
        double y = -0.176;
        int steps;
        auto f = [](double x, double y) {
            return  2 * x + tan(x + y) + (x + y + 1) * (1 + tan(x + y) * tan(x + y));
        };
        
        auto fx = [](double x, double y) {
            return  4 + 2 * tan(x + y) * tan(x + y) + 2 * (x + y + 1) * tan(x + y) * (1 + tan(x + y) * tan(x + y));
        };

        auto fy = [](double x, double y) {
            return  2 + 2 * tan(x + y) * tan(x + y) + 2 * (x + y + 1) * tan(x + y) * ( 1 +  tan(x + y) * tan(x + y));
        };


        auto g = [](double x, double y) {
            return  2 * y + tan(x + y) + (x + y + 1) * (1 + tan(x + y) * tan(x + y));
        };

        auto gx = [](double x, double y) {
            return 2 + 2 * tan(x + y) * tan(x + y) + 2 * (x + y + 1) * tan(x + y) * (1 + tan(x + y) * tan(x + y));
        };

        auto gy = [](double x, double y) {
            return 4 + 2 * tan(x + y) * tan(x + y) + 2 * (x + y + 1) * tan(x + y) * (1 + tan(x + y) * tan(x + y));
        };

        std::tie(x, y, steps) = murlib::newton_iter_exact(x, y, f, g, fx, fy, gx, gy);
        std::cout << "Min newton exact x: " << x << " y: " << y << " steps: " << steps << std::endl;

        x = -0.176;
        y = -0.176;

        std::cout << std::endl;

        std::tie(x, y, steps) = murlib::newton_iter_approx(x, y, f, g);
        std::cout << "Min newton approx x: " << x << " y: " << y << " steps: " << steps << std::endl;
    }
    return 0;
}