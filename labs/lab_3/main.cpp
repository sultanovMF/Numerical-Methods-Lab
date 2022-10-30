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
const double start_2 = 0;
const double end_2 = 2;

void print_matrix(const int n, const double* matrix) {
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
            std::cout << std::fixed << std::setprecision(2) << std::setw(6) << matrix[j * n + k] << ",\t";
        }
        std::cout << "\n";
    }
}

void print_vector(const int n, const double* vec) {
    for (int k = 0; k < n; ++k) {
        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << vec[ k] << ", ";
    }
    std::cout << "\n";
}

void build_grid_xi(double start, double end, double* test_xi, int grid_size) {
    for (int i = 0; i < grid_size; ++i) {
        test_xi[i] = start + (end - start) / (grid_size + 1) * i;
    }
}

int main() {
    // Создание тестовой сетки
    constexpr int grid_size_1 = 1000;
    const double start_1 = -3;
    const double end_1 = 3;

    double test_1_xi[grid_size_1];
    build_grid_xi(start_1, end_1, test_1_xi, grid_size_1);
    // Рациональная интерполяция
    double fracrat_test_p1_yi[grid_size_1];
    {
        int p = 1;
        int n = individual_pnum;
        double* A = new double[n * n];
        double* b = new double[n];
        murlib::build_fracrat_system(n, individual_x, individual_y, p, A, b);
        std::cout << "\n---------- Rational interpolation p = 1, q = 4 ------" << std::endl;
        print_matrix(n, A);
        std::cout << "=========== Rational interpolation val vec =========" << std::endl;
        print_vector(n, b);
        double* coef = new double[n];

        murlib::solve_gauss(n, A, b, coef);

        std::cout << "===== Rational interpolation coeff (b[0] = 1) ======" << std::endl;
        print_vector(n, coef);

        // Testing 
        for (int i = 0; i < grid_size_1; ++i) {
            fracrat_test_p1_yi[i] = murlib::fracrat_interp(n, p, coef, test_1_xi[i]);
        }
    }

    double fracrat_test_p2_yi[grid_size_1];
    {
        int p = 2;
        int n = individual_pnum;
        double* A = new double[n * n];
        double* b = new double[n];
        murlib::build_fracrat_system(n, individual_x, individual_y, p, A, b);
        std::cout << "\n---------- Rational interpolation p = 2, q = 3 ------" << std::endl;
        print_matrix(n, A);
        std::cout << "Rational interpolation val vec" << std::endl;
        print_vector(n, b);
        double* coef = new double[n];
        murlib::solve_gauss(n, A, b, coef);
        std::cout << "Rational interpolation coeff (b[0] = 1)" << std::endl;
        print_vector(n, coef);

        // Testing 
        for (int i = 0; i < grid_size_1; ++i) {
            fracrat_test_p2_yi[i] = murlib::fracrat_interp(n, p, coef, test_1_xi[i]);
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
   
    double test_grid_2[grid_size_1];
    double least_sqr_test_yi[12][grid_size_1];
    double least_sqr_point_xi[20];
    double least_sqr_point_yi[20];
    double best_r = 1000000;
    int best_deg = 1;
    int least_sqr_deg = 1;
    // Метод наименьших квадратов
    {
        build_grid_xi(0, 2, test_grid_2, grid_size_1);
        build_grid_xi(0, 2, least_sqr_point_xi, 20);
        for (int i = 0; i < 20; ++i) {
            least_sqr_point_yi[i] = individual_func(least_sqr_point_xi[i]);
        }
        for (least_sqr_deg; least_sqr_deg <= 12; ++least_sqr_deg) {
            double* A = new double[least_sqr_deg * least_sqr_deg];
            double* b = new double[least_sqr_deg];
            double* S = new double[least_sqr_deg * least_sqr_deg];
            double* coef = new double[least_sqr_deg];
            murlib::build_least_sqr_sys(20, least_sqr_deg, least_sqr_point_xi, least_sqr_point_yi, A, b);
           // std::cout << "\n======= Least squares symmetric matrix =======" << std::endl;
           // print_matrix(least_sqr_deg, A);
           // std::cout << std::endl;
            //print_vector(least_sqr_deg, b);
            murlib::sqroot_decompostion(least_sqr_deg, A, S);
            murlib::sqroot_solve(least_sqr_deg, S, b, coef);
           // std::cout << "\n======= Least squares coeffs =======" << std::endl;
            //print_vector(least_sqr_deg, coef);

            // Testing 
            for (int i = 0; i < grid_size_1; ++i) {
                least_sqr_test_yi[least_sqr_deg - 1][i] = murlib::least_sqr_approx(least_sqr_deg, coef, test_grid_2[i]);
            }
            // Find max error
            double current_max = 0;
            for (int i = 0; i < 20; ++i) {
                double R = std::abs(murlib::least_sqr_approx(least_sqr_deg, coef, least_sqr_point_xi[i]) - least_sqr_point_yi[i]);
                if (R > current_max) {
                    current_max = R;
                }
            }
            if (current_max < best_r) {
                best_r = current_max;
                best_deg = least_sqr_deg;
            }
            delete[] A;
            delete[] b;
            delete[] S;
            delete[] coef;
        }
    }
    least_sqr_deg = 5;
    std::cout << "\n======= Least squares =======" << std::endl;
    std::cout << "Best degree = " << best_deg << " with error " << std::setprecision(12) << best_r << std::endl;
    std::cout << std::endl;
    
    // Пятидиагональная прогонка
    // Скорее всего неверно составлена система или выводится матрица
    using Random = effolkronium::random_static;
    {
        const int n = 10;
        double M[n * n];
        std::fill(M, M + n * n, 0);
        double A[n];
        double B[n];
        double C[n];
        double D[n];
        double E[n];
        double F[n];

        A[0] = 0;
        A[1] = 0;
        A[n-1] = Random::get(-10, 10);
        A[n-2] = Random::get(-10, 10);

        B[0] = 0;
        B[1] =   Random::get(-10, 10);
        B[n-1] = Random::get(-10, 10);
        B[n-2] = Random::get(-10, 10);

        C[0] =   Random::get(-10, 10);
        C[1] =   Random::get(-10, 10);
        C[n-1] = Random::get(-10, 10);
        C[n-2] = Random::get(-10, 10);

        D[0] =   Random::get(-10, 10);
        D[1] =   Random::get(-10, 10);
        D[n - 1] = 0;
        D[n-2] = Random::get(-10, 10);

        E[0] = Random::get(-10, 10);
        E[1] = Random::get(-10, 10);
        E[n - 1] = 0;
        E[n - 2] = 0;

        F[0] = Random::get(-10, 10);
        F[1] = Random::get(-10, 10);
        F[n-1] = Random::get(-10, 10);
        F[n-2] = Random::get(-10, 10);
        for (int i = 2; i < n-2; ++i) {
            A[i] = Random::get(-10, 10);
            B[i] = Random::get(-10, 10);
            C[i] = Random::get(-10, 10);
            D[i] = Random::get(-10, 10);
            E[i] = Random::get(-10, 10);
            F[i] = Random::get(-10, 10);
        }
        double X[n];

        murlib::pentadiagonal_matrix(n, A, B, C, D, E, F, X);

        std::cout << "\n======= Pentadiagonal matrix =======" << std::endl;
        M[0] = C[0];
        M[1] = D[0];
        M[2] = E[0];

        M[n] = B[1];
        M[n + 1] = C[1];
        M[n + 2] = D[1];
        M[n + 3] = E[1];

        for (int i = 2; i < n-2; ++i) {
            M[i * n + i - 2] = A[i];
            M[i * n + i - 1] = B[i];
            M[i * n + i] = C[i];
            M[i * n + i + 1] = D[i];
            M[i * n + i + 1] = E[i];
        }

        M[n * (n - 1) - 4] = A[n-2];
        M[n * (n - 1) - 3] = B[n-2];
        M[n * (n - 1) - 2] = C[n-2];
        M[n * (n - 1) - 1] = D[n-2];

        M[n * n - 3] = A[n - 1];
        M[n * n - 2] = B[n - 1];
        M[n * n - 1] = C[n - 1];

        print_matrix(n, M);
        std::cout << "RHS: ";
        print_vector(n, F);
        std::cout << "SOL: ";
        print_vector(n, X);
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


        ImGui::Begin("Lab 3");

        
        if (ImPlot::BeginPlot("p = 1, q = 4")) {
            ImPlot::PlotLine("I(x)", test_1_xi, fracrat_test_p1_yi, grid_size_1, ImPlotLineFlags_SkipNaN);
            ImPlot::PlotScatter("Interpolation points", individual_x, individual_y, individual_pnum);
            ImPlot::EndPlot();
        }

        if (ImPlot::BeginPlot("p = 2, q = 3")) {
            ImPlot::PlotLine("I(x)", test_1_xi, fracrat_test_p2_yi, grid_size_1, ImPlotLineFlags_SkipNaN);
            ImPlot::PlotScatter("Interpolation points", individual_x, individual_y, individual_pnum);
            ImPlot::EndPlot();
        }
        ImGui::SliderInt("Degree of polynom", &least_sqr_deg, 1, 12);
        if (ImPlot::BeginPlot("Least square")) {
            ImPlot::PlotLine("f(x)", test_grid_2, least_sqr_test_yi[least_sqr_deg-1], grid_size_1, ImPlotLineFlags_SkipNaN);
            ImPlot::PlotScatter("Interpolation points", least_sqr_point_xi, least_sqr_point_yi, 20);
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