#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>

// Самописная библиотека мурлиб для классических численных методов
#include <murlib/interpolation.h>
#include <murlib/matrix_solver.h>
#include <murlib/approximation.h>
#include <murlib/miscellaneous.h>

// Библиотеки для вывода графиков/табличек
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <imgui-SFML.h>
#include <imgui.h>
#include <implot/implot.h>

// Библиотека для рандома
#include <effolkronium/random.hpp>

//* Индивидуальное задание
//  Колличество точек интерполяции
constexpr int individual_pnum = 6;
// Табличные значения аргумента
constexpr double individual_x[] = {
	-2.25, -1.65, -1.05, -0.45, 0.15, 0.75
};
// Табличные значения функции
constexpr double individual_y[] = {
	0.29422, 0.18737, -1.0215, -5.4471, -1.4440, 2.5873
};
// Функция из первой лабораторной
double individual_func(double x) {
	return std::atan(x) / (1 + x * x);
};
constexpr int grid_size = 1000;
const double start = -3;
const double end = 3;

void print_matrix(const int n, const double* matrix) {
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            std::cout << std::fixed << std::setprecision(6) << matrix[j * n + k] << ",\t";
        }
        std::cout << "\n";
    }
}

void print_vector(const int n, const double* vec) {
    for (int k = 0; k < n; ++k) {
        std::cout << std::fixed << std::setprecision(6) << vec[ k] << ", ";
    }
    std::cout << "\n";
}
int main() {
    // Создание тестовой сетки
    double test_xi[grid_size];
    for (int i = 0; i < grid_size; ++i) {
        test_xi[i] = start + (end - start) / (grid_size + 1)*i;
    }
    // Рациональная интерполяция
    std::cout << "----------- p = 1, q = 4 -------------" << std::endl;
    double fracrat_test_p1_yi[grid_size];
    {
        int p = 1;
        int n = individual_pnum;
        double* A = new double[n * n];
        double* b = new double[n];
        murlib::build_fracrat_system(n, individual_x, individual_y, p, A, b);
        std::cout << "=========== Rational interpolation matrix ==========" << std::endl;
        print_matrix(n, A);
        std::cout << "=========== Rational interpolation val vec =========" << std::endl;
        print_vector(n, b);
        double* coef = new double[n];

        murlib::solve_gauss(n, A, b, coef);

        std::cout << "===== Rational interpolation coeff (b[0] = 1) ======" << std::endl;
        print_vector(n, coef);

        // Testing 
        for (int i = 0; i < grid_size; ++i) {
            fracrat_test_p1_yi[i] = murlib::fracrat_interp(n, p, coef, test_xi[i]);
        }
    }
    std::cout << "\n----------- p = 2, q = 3 -------------" << std::endl;
    double fracrat_test_p2_yi[grid_size];
    {
        int p = 2;
        int n = individual_pnum;
        double* A = new double[n * n];
        double* b = new double[n];
        murlib::build_fracrat_system(n, individual_x, individual_y, p, A, b);
        std::cout << "=========== Rational interpolation matrix ==========" << std::endl;
        print_matrix(n, A);
        std::cout << "=========== Rational interpolation val vec =========" << std::endl;
        print_vector(n, b);
        double* coef = new double[n];
        murlib::solve_gauss(n, A, b, coef);
        std::cout << "===== Rational interpolation coeff (b[0] = 1) ======" << std::endl;
        print_vector(n, coef);

        // Testing 
        for (int i = 0; i < grid_size; ++i) {
            fracrat_test_p2_yi[i] = murlib::fracrat_interp(n, p, coef, test_xi[i]);
        }
    }
    std::cout << "\n----------- LU p = 2, q = 3 -------------" << std::endl;
    //double fracrat_test_p2_lu_yi[grid_size];
    {
        int p = 2;
        int n = individual_pnum;
        double* A = new double[n * n];
        double* b = new double[n];
        murlib::build_fracrat_system(n, individual_x, individual_y, p, A, b);
        double* coef = new double[n];

        double* L = new double[n * n];
        double* U = new double[n * n];

        murlib::lu_decompostion(n, A, L, U);
        murlib::solve_lu(n, L, U, b, coef);
        std::cout << "===== Rational interpolation coeff (b[0] = 1) ======" << std::endl;
        print_vector(n, coef);

        // Testing 
        //for (int i = 0; i < grid_size; ++i) {
        //    fracrat_test_p2_lu_yi[i] = murlib::fracrat_interp(n, p, coef, test_xi[i]);
        //}

        delete[] L;
        delete[] U;
    }
   
    // Метод наименьших квадратов


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


        ImGui::Begin("Lab 3");

        
        if (ImPlot::BeginPlot("p = 1, q = 4")) {
            ImPlot::PlotLine("I(x)", test_xi, fracrat_test_p1_yi, grid_size, ImPlotLineFlags_SkipNaN);
            ImPlot::PlotScatter("Interpolation points", individual_x, individual_y, individual_pnum);
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("p = 2, q = 3")) {
            ImPlot::PlotLine("I(x)", test_xi, fracrat_test_p2_yi, grid_size, ImPlotLineFlags_SkipNaN);
            ImPlot::PlotScatter("Interpolation points", individual_x, individual_y, individual_pnum);
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