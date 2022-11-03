#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>
#include <array>

// Самописная библиотека мурлиб для классических численных методов
#include <murlib/integration.h>
#include <murlib/miscellaneous.h>

// Библиотеки для вывода графиков/табличек
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <implot/implot.h>



const double true_value = std::atan(2) * std::atan(2) / 2;
constexpr double delta = 1e-3;

double individual_func(double x) {
	return std::atan(x) / (1 + x * x);
};

double t_to_x(double start, double end, double t) {
    return (start + end) / 2 + (end - start) / 2 * t;
}

double individual_func_t(double start, double end, double t) {
    double x = t_to_x(start, end, t);
    return std::atan(x) / (1 + x * x);
}

void build_grid(double start, double end, double* test_xi, int grid_size) {
	for (int i = 0; i < grid_size; ++i) {
		test_xi[i] = start + (end - start) / (grid_size + 1) * i;
	}
}


constexpr double a = 0;
constexpr double b = 2;

int main() {
    // Классические формулы
	double rectangle_error[99];
	double trapezoid_error[99];
	double simpson_error  [99];
    std::cout << "True value: " << true_value << std::endl;
    int rectangle_optimal_n = 0;
    int trapezoid_optimal_n = 0;
    int simpson_optimal_n = 0;
    {

        int rectangle_optimal_i = 0;
        int trapezoid_optimal_i = 0;
        int simpson_optimal_i = 0;

        double rectangle_optimal_value = 0;
        double trapezoid_optimal_value = 0;
        double simpson_optimal_value = 0;
        int i = 0;
        for (int n = 100; n < 10e3; n += 100, ++i) {
            
            double* test_xi = new double[n];
            build_grid(a, b, test_xi, n);
            double rectangle_value = murlib::integrate_rectangle(test_xi, individual_func, n);
            double trapezoid_value = murlib::integrate_trapezoid(test_xi, individual_func, n);
            double simpson_value = murlib::integrate_simpson(test_xi, individual_func, n);

            rectangle_error[i] = std::abs(rectangle_value - true_value);
            trapezoid_error[i] = std::abs(trapezoid_value - true_value);
            simpson_error[i] = std::abs(simpson_value - true_value);

            if ((rectangle_error[i] < delta) && (rectangle_optimal_n == 0)) {
                rectangle_optimal_n = n;
                rectangle_optimal_i = i;
                rectangle_optimal_value = rectangle_value;
            }

            if ((trapezoid_error[i] < delta) && (trapezoid_optimal_n == 0)) {
                trapezoid_optimal_n = n;
                trapezoid_optimal_i = i;
                trapezoid_optimal_value = trapezoid_value;
            }
            if ((simpson_error[i]   < delta) && (simpson_optimal_n == 0)) {
                simpson_optimal_n = n;
                simpson_optimal_i = i;
                simpson_optimal_value = simpson_value;
            }
           //std::cout << i << std::endl;
            delete[] test_xi;
        }
      //  0.001
        std::cout << i << std::endl;
        std::cout << "Optimal rectangle n = " << rectangle_optimal_n
            << " with error: " << rectangle_error[rectangle_optimal_i]
            << "\tand value: " << rectangle_optimal_value
            <<  std::endl;
        std::cout << "Optimal trapezoid n = " << trapezoid_optimal_n
            << " with error: " << trapezoid_error[trapezoid_optimal_i]
            << "\tand value: " << trapezoid_optimal_value
            << std::endl;
        std::cout << "Optimal simpson   n = " << simpson_optimal_n 
            << "  with error: " << simpson_error[simpson_optimal_i] 
            << "\tand value: " << simpson_optimal_value
            << std::endl;
    }

    std::cout << std::endl;
    // Процесс Эйткена

    {
        const int n = 10;
        double p = murlib::aitken_rectangle(a, b, n, individual_func);
        std::cout << "Rectangle optimal p = " << p << std::endl;

        p = murlib::aitken_trapezoid(a, b, n, individual_func);
        std::cout << "Trapezoid optimal p = " << p << std::endl;

        p = murlib::aitken_simpson(a, b, n, individual_func);
        std::cout << "Simpson optimal p = " << p << std::endl;
    }

    // Метод Рунге
    int p = round(murlib::aitken_rectangle(0, 2, rectangle_optimal_n, individual_func));
    {
        double runge_rectangle_res = murlib::runge_rectangle(0, 2, rectangle_optimal_n, p, individual_func);
        std::cout << "Rectangle runge: " << runge_rectangle_res << std::endl;
        std::cout << "Rectangle runge error: " << std::fixed << std::setprecision(12) << abs(runge_rectangle_res - true_value) << std::endl;
    }

    //Метод Ромберга
    {
        for (int q = 1; q < 10; ++q) {
            double romberg_rectangle_res = murlib::romberg_rectangle(0, 2, rectangle_optimal_n, p, q, individual_func);
            std::cout << "Romberg with q = " << q << " " << std::setprecision(20) << romberg_rectangle_res << " and error "  << abs(romberg_rectangle_res - true_value)  << std::endl;
        }
    }
    // Квадратурные формулы Гаусса TODO
    double gauss_error[6];
    {
        double gauss[6];
        
        gauss[0] = 2.0 * individual_func_t(a, b, 0.);
       
        gauss[1] = 1.0 * individual_func_t(a, b, -0.5773502692) 
                + 1.0 * individual_func_t(a, b, 0.5773502692);
        
        gauss[2] = 0.5555555556 * individual_func_t(a, b, -0.7745966692) 
                + 0.8888888888 * individual_func_t(a, b,  0.)
                + 0.5555555556 * individual_func_t(a, b,  0.7745966692);
        
        gauss[3] = 0.3478548451 * individual_func_t(a, b, -0.8611363115)
                + 0.6521451549 * individual_func_t(a, b, -0.3399810436)
                + 0.6521451549 * individual_func_t(a, b,  0.3399810436)
                + 0.3478548451 * individual_func_t(a, b,  0.8611363115);
        
        gauss[4] = 0.2369268851 * individual_func_t(a, b, -0.9061798459)
                + 0.4786286705 * individual_func_t(a, b, -0.5384693101)
                + 0.5688888888 * individual_func_t(a, b,  0.0)
                + 0.4786286705 * individual_func_t(a, b,  0.5384693101)
                + 0.2369268851 * individual_func_t(a, b,  0.9061798459);
        
        gauss[5] = 0.1713244924 * individual_func_t(a, b, -0.9324695142)
                + 0.3607615730 * individual_func_t(a, b, -0.6612093864)
                + 0.4679139346 * individual_func_t(a, b, -0.2386191861)
                + 0.4679139346 * individual_func_t(a, b,  0.2386191861)
                + 0.3607615730 * individual_func_t(a, b,  0.6612093864)
                + 0.1713244924 * individual_func_t(a, b,  0.9324695142);
        
        for (int i = 0; i < 6; ++i) {
            gauss[i] *= (b - a) / 2.;
            gauss_error[i] = abs(gauss[i] - true_value);
            std::cout << "Gauss n = " << i + 1 << " with error " << gauss_error[i] << std::endl;
            //if (gauss_error[i] < delta) {
            //    std::cout << "Gauss optimal" << i + 1 << std::endl;
            //} 4
        }
    }
    // Адаптивная сетка
    double* adaptive_grid_xi;
    double adaptive_value = 0;
    int adaptive_n = 2;
    {
        // murlib::adaptive_integration();
        for (;;adaptive_n++) {
            double* xi = new double[adaptive_n];
            adaptive_value = murlib::adaptive_integration(delta, a, b, adaptive_n, xi, individual_func);
            if (abs(adaptive_value - true_value) < delta) {
                adaptive_grid_xi = new double[adaptive_n];
                std::copy(xi, xi + adaptive_n, adaptive_grid_xi);
                break;
            }
            delete[] xi;
        }
    }
    adaptive_n--;

    std::cout << std::endl << "Adaptive n: " << adaptive_n <<" value " << adaptive_value << " with error " << abs(adaptive_value - true_value) << " (" << delta << ")" << std::endl;
    double* adaptive_grid_yi = new double[adaptive_n];
    for (int i = 0; i < adaptive_n; ++i) {
        adaptive_grid_yi[i] = individual_func(adaptive_grid_xi[i]);
    }

    // Монте-Карло
    std::cout << "Monte Carlo" << std::endl;
    double monte_carlo_estimate[97];
    {
        for (int n = 2; n < 100; ++n) {
            double result = murlib::monte_carlo(a, b, n, individual_func);
          //  std::cout << result << std::endl;
            monte_carlo_estimate[n - 2] = result / (b - a);
        }
    }
    sf::RenderWindow window(sf::VideoMode(1280, 720), "ImGui + SFML = <3");
    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);

    sf::Clock deltaClock;

    ImGui::CreateContext();
    ImPlot::CreateContext();


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(event);

            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());

        static bool use_work_area = true;
        static ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
        // ImGui::ShowDemoWindow();
        const ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(use_work_area ? viewport->WorkPos : viewport->Pos);
        ImGui::SetNextWindowSize(use_work_area ? viewport->WorkSize : viewport->Size);


        ImGui::Begin("Lab 2");
        if (ImPlot::BeginPlot("Rectangle error")) {
            ImPlot::PlotBars("R", rectangle_error, 99);
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Trapezoid error")) {
            ImPlot::PlotBars("T", trapezoid_error, 99);
            ImPlot::EndPlot();
        }
        if (ImPlot::BeginPlot("Simpson error")) {
            ImPlot::PlotBars("S", simpson_error, 99);
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Adaptive distribution")) {
            ImPlot::PlotScatter("Adaptive points", adaptive_grid_xi, adaptive_grid_yi, adaptive_n);
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Monte-Carlo estimate")) {
            ImPlot::PlotBars("Estimate", monte_carlo_estimate, 97);
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("Gauss estimate")) {
            ImPlot::PlotBars("Estimate", gauss_error, 6);
            ImPlot::EndPlot();
        }

        ImGui::End();

        window.clear();
        ImGui::SFML::Render(window);
        window.display();
    }


    ImGui::SFML::Shutdown();

    delete[] adaptive_grid_xi;
    delete[] adaptive_grid_yi;
	return 0;
}