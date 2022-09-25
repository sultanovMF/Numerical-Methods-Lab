#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <iomanip>
#include <chrono>

#include <murlib/interpolation.h>

#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <implot/implot.h> 

const double a = 0;
const double b = 2;
const double delta = 10e-3;
double individual_func(double x) {
	return atan(x) / (1 + x * x);
};

const int grid_size = 10e2;


class Polynom {
public:
	Polynom() {}
	Polynom(int degree, std::vector<double>  grid_x, std::vector<double>  grid_y) : degree(degree), grid_x(grid_x), grid_y(grid_y) {}

	double lagrange_interpolation(double t) {
		return murlib::lagrange_polynom(grid_x, grid_y, t, degree);
	}

    double newton_interpolation(double t) {
        return murlib::newton_polynom(grid_x, grid_y, t, degree);
    }

private:
	unsigned short degree;

	std::vector<double> grid_x;
	std::vector<double> grid_y;

};

void build_uniform_grid(std::vector<double>& grid_x, std::vector<double>& grid_y, const int degree, const double start, const double end) {
    for (int i = 0; i < degree; ++i) {
        // Filling data into interpolation nodes
        grid_x[i] = start + (end - start) / degree * i;
        grid_y[i] = individual_func(grid_x[i]);
    }
}

int main() {
	std::array <Polynom, 15> lagrange;
	for (int n = 1; n <= 15; ++n) {
		std::vector<double> grid_x(n);
		std::vector<double> grid_y(n);

        build_uniform_grid(grid_x, grid_y, n, a, b);

		lagrange[n - 1] = Polynom(n, grid_x, grid_y);
	}


	std::array <double, 15> max_error;
	std::array <double, grid_size> x;
	std::array <double, grid_size> f;
	std::array <std::array<double, grid_size>, 15> L;
	std::array <std::array<double, grid_size>, 15> error;

	for (int i = 0; i < grid_size; ++i) {
		x[i] = a + (b - a) / grid_size * i;
		f[i] = individual_func(x[i]);
		for (int n = 1; n <= 15; ++n) {
			L[n-1][i] = lagrange[n-1].lagrange_interpolation(x[i]);
			error[n - 1][i] = abs(L[n - 1][i] - f[i]);
		}
	}

	for (short n = 1; n <= 15; ++n) {
		max_error[n-1] = *std::max_element(error[n-1].begin(), error[n-1].end());
	}
	auto min_error = std::min_element(max_error.begin(), max_error.end());
	int min_error_degree = std::distance(max_error.begin(), min_error) + 1;
	double min_error_val = *min_error;

    std::vector<double> newton_grid_x(min_error_degree);
    std::vector<double> newton_grid_y(min_error_degree);
    build_uniform_grid(newton_grid_x, newton_grid_y, min_error_degree, a, b);

    Polynom weird_polynom(min_error_degree, newton_grid_x, newton_grid_y);

    std::array <double, grid_size> N_n0;
    double newton_max_error = 0;

    for (int i = 0; i < grid_size; ++i) {
        N_n0[i] = weird_polynom.newton_interpolation(a + (b - a) / grid_size * i);
        if (abs(N_n0[i] - f[i]) > newton_max_error) {
            newton_max_error = abs(N_n0[i] - f[i]);
        }
    }
    int degree_of_polynom = min_error_degree;
    sf::RenderWindow window(sf::VideoMode(1280, 720), "ImGui + SFML = <3");
    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);

    sf::Clock deltaClock;

    ImGui::CreateContext();
    ImPlot::CreateContext();

    bool show_newton_polynom;


   // /* Speed Test */
   // auto begin = std::chrono::steady_clock::now();
   // for (int i = 0; i < grid_size; i++) {
   //     weird_polynom.lagrange_interpolation(a + (b - a) / grid_size * i);
   // }
   // auto end = std::chrono::steady_clock::now();
   // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
   // begin = std::chrono::steady_clock::now();
   // for (int i = 0; i < grid_size; i++) {
   //     weird_polynom.newton_interpolation(a + (b - a) / grid_size * 2);
   // }
   // end = std::chrono::steady_clock::now();
   // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
   //// auto lagrange_time = duration_cast<std::chrono::milliseconds>(stop - start);


    //start = std::chrono::high_resolution_clock::now();

    //stop = std::chrono::high_resolution_clock::now();
  //  auto newton_time = duration_cast<std::chrono::milliseconds>(stop - start);

    //std::cout << newton_time << " " << lagrange_time << std::endl;
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(event);

            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        ImGui::SFML::Update(window, deltaClock.restart());
        
        ImGui::Begin("My Window");
        ImGui::Checkbox("show newton", &show_newton_polynom);
        ImGui::SliderInt("degree of polynom", &degree_of_polynom, 1, 15);
        if (show_newton_polynom) {
          //  ImGui::Text("Newton speed: %lld ns \t Lagrange speed: %lld ns", newton_time, lagrange_time);
        }
        if (ImPlot::BeginPlot("Function and its lagrange interpolation")) {
            ImPlot::PlotLine("f(x)", x.data(), f.data(), grid_size);
            ImPlot::PlotLine("L_n(x)", x.data(), L[degree_of_polynom-1].data(), grid_size);
           // ImPlot::PlotLine("N_n(x)", x.data(), N.data(), grid_size);
            if (show_newton_polynom) {
                ImPlot::PlotLine("N_n(x)", x.data(), N_n0.data(), grid_size);
            }
            ImPlot::EndPlot();
        }


        if (ImGui::BeginTable("error", 2))
        {
            ImGui::TableNextColumn();

            if (ImGui::BeginTable("error_table", 2))
            {
                ImGui::TableNextColumn();
                ImGui::Text("n");
                ImGui::TableNextColumn();
                ImGui::Text("error");
                for (int n = 1; n <= 15; ++n) {
                    ImGui::TableNextColumn();
                    ImGui::Text("%i", n);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.15f", max_error[n - 1]);
                }
                ImGui::TableNextColumn();
                ImGui::Text(("Newton error at " + std::to_string(min_error_degree)).c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.15f", newton_max_error);
                ImGui::EndTable();
            }

            ImGui::TableNextColumn();

            if (ImPlot::BeginPlot("max_error")) {
                ImPlot::PlotBars("max_error", max_error.data(), 15);
                ImPlot::EndPlot();
            }


            ImGui::EndTable();


        }

        if (ImPlot::BeginPlot("Error(n_0) function")) {
            ImPlot::PlotLine("f(x)", x.data(), error[min_error_degree-1].data(), grid_size);
            ImPlot::EndPlot();
        }

        ImGui::End();
        
        window.clear();
        ImGui::SFML::Render(window);
        window.display();
    }


    ImGui::SFML::Shutdown();
    return 0;
}


