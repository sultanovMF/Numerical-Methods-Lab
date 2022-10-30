#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>
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




//* Константы и функции индивидуального задания:
//* Отрезок, на котором происходит интерполяция/аппроксимация
const double a = 0, b = 2;
//* Рекомендованное значение ошибки
const double delta = 10e-3;
//* Размер контрольной сетки
const int grid_size = 1000;
//* Максимальная степень многочлен
const int max_points_number = 15;
//* Функция, которую требуется аппроксимировать
double individual_func(double x) {
    return std::atan(x) / (1 + x * x);
};

double normilize(const double x) {
    return 2 * murlib::PI / (b - a) * (x - a);
}

double unnormilize(const double x) {
    return (b - a) / (2 * murlib::PI) * x + a;
}

double P(double x) {
    // deg P = 20
    return 0.5353981636e0 - 0.1426990818 * x - 0.1786504591 * pow
    (x - 0.1e1, 2) + 0.2916666667e0 * pow(x - 0.1e1, 3) -
        0.2023414371e0 * pow(x - 0.1e1, 4) + 0.4400810377e-1 * pow(x
            - 0.1e1, 5) + 0.6757928148e-1 * pow(x - 0.1e1, 6) -
        0.9404761905e-1 * pow(x - 0.1e1, 7) + 0.6025797832e-1 * pow
        (x - 0.1e1, 8) - 0.1149805769e-1 * pow(x - 0.1e1, 9) -
        0.2019343148e-1 * pow(x - 0.1e1, 10) + 0.2665268759e-1 * pow
        (x - 0.1e1, 11) - 0.1655597185e-1 * pow(x - 0.1e1, 12) +
        0.2929147289e-2 * pow(x - 0.1e1, 13) + 0.5627856495e-2 * pow
        (x - 0.1e1, 14) - 0.7222638473e-2 * pow(x - 0.1e1, 15) +
        0.4408710225e-2 * pow(x - 0.1e1, 15) - 0.7399461360e-3 * pow
        (x - 0.1e1, 17) - 0.1518662449e-2 * pow(x - 0.1e1, 18) +
        0.1914334530e-2 * pow(x - 0.1e1, 19);
}

int main() {

    // Создание тестовой сетки
    double test_xi[grid_size];
    double target_yi[grid_size];
    double P_yi[grid_size];
    for (int i = 0; i < grid_size; ++i) {
        test_xi[i] = a + (b - a) / (grid_size) * i;
        target_yi[i] = individual_func(test_xi[i]);
        P_yi[i] = P(test_xi[i]);
    }
    // МНРП!
    double mnrp_result_yi[grid_size];
    double mnrp_error[grid_size];
    int mnrp_best_degree = 20;
    for (int i = mnrp_best_degree; i >= 0; --i) {
        murlib::mnrp_for_polynoms(P, i, a, b, grid_size, test_xi, mnrp_result_yi, mnrp_error);
        for (int i = 0; i < grid_size; ++i) {
            P_yi[i] = mnrp_result_yi[i];
        }
        double mnrp_el= *std::max_element(mnrp_error, mnrp_error + grid_size);


        if (mnrp_el > delta) {
            mnrp_best_degree = i + 1;
            break;
        }
    }
    // Интерполяция сплайнами
    int spline_count = 2;
    double spline_test_yi[grid_size];
    double spline_error = 10000;
    double spline_error_data[grid_size];
    while (spline_error > delta && spline_count < 100) {
        double* A = new double[spline_count];
        double* B = new double[spline_count];
        double* C = new double[spline_count];
        double* D = new double[spline_count];

        double* spline_xi = new double[spline_count + 1];
        double* spline_yi = new double[spline_count + 1];

        for (int i = 0; i <= spline_count; ++i) {
            spline_xi[i] = a + (b - a) / (spline_count)*i;
            spline_yi[i] = individual_func(spline_xi[i]);
        }
        murlib::build_splines(spline_count, A, B, C, D, spline_xi, spline_yi);

        //std::cout << A[0] << " " << B[0] << " " << C[0] << " " << D[0] << std::endl;
        //std::cout << A[1] << " " << B[1] << " " << C[1] << " " << D[1] << std::endl;

        int spline_zone = 0;
        double spline_current_max_error = 0;
        for (int i = 0; i < grid_size; ++i) {
            spline_test_yi[i] = murlib::spline_iterpolation(spline_count, A, B, C, D, spline_xi, spline_zone, test_xi[i]);
            spline_error_data[i] = std::abs(spline_test_yi[i] - target_yi[i]);
            if (std::abs(spline_test_yi[i] - target_yi[i]) > spline_current_max_error) {
                spline_current_max_error = std::abs(spline_test_yi[i] - target_yi[i]);
                
            }

            if (test_xi[i + 1] > (a + (b-a) / spline_count * (spline_zone + 1)) && (i < grid_size - 1)) {
                spline_zone++;
            }

        } 
        spline_error = spline_current_max_error;
      //  std::cout << spline_error << std::endl;
        delete[] A;
        delete[] B;
        delete[] C;
        delete[] D;
        delete[] spline_xi;
        delete[] spline_yi;
        spline_count++;
    }


    /// Интерполяция методом Лагранжа
    double lagrange_xi[max_points_number][max_points_number];
    double lagrange_yi[max_points_number][max_points_number];

    for (int n = 1; n <= max_points_number; ++n) {
        for (int i = 0; i < n; ++i) {
            lagrange_xi[n - 1][i] = a + (b - a) / (n)*i;
            lagrange_yi[n - 1][i] = individual_func(lagrange_xi[n - 1][i]);
        }
    }

    double lagrange_test_yi[max_points_number][grid_size];

    double lagrange_error[max_points_number][grid_size];
    int lagrange_best_error_index;

    for (int i = 0; i < grid_size; ++i) {
        for (int n = 1; n <= max_points_number; ++n) {
            lagrange_test_yi[n-1][i] = murlib::lagrange_polynom(lagrange_xi[n-1], lagrange_yi[n - 1], test_xi[i], n); //n+1
            lagrange_error[n-1][i] = std::abs(lagrange_test_yi[n-1][i] - target_yi[i]);
        }
    }

    double lagrange_max_error[max_points_number];
    std::cout << "Lagrange error" << std::endl;
    for (int n = 0; n < max_points_number; n++) {
        double* max_element = std::max_element(lagrange_error[n], lagrange_error[n] + grid_size);
        lagrange_max_error[n] = *max_element;
        std::cout << n+1 << "\t" << lagrange_max_error[n] << std::endl;
    }
   
    {
        // Нахождение индекса максимального элемента
        auto min_element = std::min_element(lagrange_max_error, lagrange_max_error + max_points_number);
        lagrange_best_error_index =  std::distance(lagrange_max_error, min_element);
    }

    double* chebyshev_grid_xi = new double[lagrange_best_error_index + 1];
    double* chebyshev_grid_yi = new double[lagrange_best_error_index + 1];

    for (int i = 0; i < lagrange_best_error_index + 1; ++i) {
        chebyshev_grid_xi[i] = murlib::get_chebyshev_root(a, b, lagrange_best_error_index + 1, i);
        chebyshev_grid_yi[i] = individual_func(chebyshev_grid_xi[i]);
    }
    double chebyshev_test_yi[grid_size];
    double chebyshev_max_error = 0;
    for (int i = 0; i < grid_size; ++i) {
        chebyshev_test_yi[i] = murlib::lagrange_polynom(chebyshev_grid_xi, chebyshev_grid_yi, test_xi[i], lagrange_best_error_index + 1);
        if (std::abs(chebyshev_test_yi[i] - target_yi[i]) > chebyshev_max_error) {
            chebyshev_max_error = std::abs(chebyshev_test_yi[i] - target_yi[i]);
        }
    }

    ///// Интерполяция полиномом Ньютона (разделенные разности)
    double newton_test_yi[grid_size];
    double newton_error = 0;
    {
        const int points_number = lagrange_best_error_index + 1;
        double** v = new double* [points_number];
        for (int i = 0; i < points_number; ++i) {
            v[i] = new double[points_number];
            v[0][i] = lagrange_yi[lagrange_best_error_index][i];
        }

        for (int i = 1; i < points_number; ++i) {
            for (int j = 0; i < points_number; ++i) {
                v[i][j] = 0;
            }
        }

        for (int i = 0; i < grid_size; ++i) {
            newton_test_yi[i] = murlib::newton_polynom(
                lagrange_xi[lagrange_best_error_index], v, test_xi[i], points_number
            );

            if (std::abs(newton_test_yi[i] - target_yi[i]) > newton_error) {
                newton_error = std::abs(newton_test_yi[i] - target_yi[i]);
            }
        }

        for (int i = 0; i < points_number; ++i) {
            delete[] v[i];
        }
        delete[] v;
        //return murlib::newton_polynom(grid_x, v, t, points_number);
    }
    
    /// Тригонометрическая интерполяция
    double trigonometry_error = 10000;
    double trigonometry_test_yi[grid_size];
    double trigonometry_error_data[grid_size];

    int trig_degree = 3;
    while (trigonometry_error > delta && trig_degree < 21) {
        const int trig_point_number = trig_degree * 2 + 1;
        double* trigonometry_xi = new double[trig_point_number];
        double* trigonometry_yi = new double[trig_point_number];

        for (int i = 0; i < trig_point_number; ++i) {
            trigonometry_xi[i] = i * (2 * murlib::PI) / (trig_point_number);
            trigonometry_yi[i] = individual_func(unnormilize(trigonometry_xi[i]));
        }
        double trigonometry_current_error = 0;
        for (int i = 0; i < grid_size; ++i) {
            trigonometry_test_yi[i] = murlib::trigonometric_polynom(normilize(test_xi[i]), trigonometry_xi, trigonometry_yi, trig_degree);
            trigonometry_error_data[i] = (std::abs(trigonometry_test_yi[i] - target_yi[i]));
            if ((std::abs(trigonometry_test_yi[i] - target_yi[i]) > trigonometry_current_error) && (test_xi[i] < unnormilize(trigonometry_xi[trig_point_number-1]))) {
                trigonometry_current_error = std::abs(trigonometry_test_yi[i] - target_yi[i]);
            }
        }

        delete[] trigonometry_xi;
        delete[] trigonometry_yi;
        trigonometry_error = trigonometry_current_error;
       // std::cout << trig_degree << " " << trigonometry_error << std::endl;
        trig_degree++;
    }

    //const int trig_degree = 201;
    //const int trig_point_number = trig_degree * 2 + 1;
    //double trigonometry_xi[trig_point_number];
    //double trigonometry_yi[trig_point_number];
    //double trigonometry_test_yi[grid_size];

    //for (int i = 0; i < trig_point_number; ++i) {
    //trigonometry_xi[i] = i * (2 * murlib::PI) / (trig_point_number);
    //// TODO спросить про i * (2 * murlib::PI) / (trig_point_number-1);
    //trigonometry_yi[i] = individual_func(unnormilize(trigonometry_xi[i]));
    //}


    //for (int i = 0; i < grid_size; ++i) {
    //    trigonometry_test_yi[i] = murlib::trigonometric_polynom(normilize(test_xi[i]), trigonometry_xi, trigonometry_yi, trig_degree);
    //}
    sf::RenderWindow window(sf::VideoMode(1280, 720), "ImGui + SFML = <3");
    window.setFramerateLimit(60);
    ImGui::SFML::Init(window);

    sf::Clock deltaClock;

    ImGui::CreateContext();
    ImPlot::CreateContext();

    int current_plot_degree = lagrange_best_error_index;

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
        

        ImGui::Begin("Lab 1");
        if (ImGui::BeginTable("error", 2))
        {
            ImGui::TableNextColumn();

            if (ImGui::BeginTable("error_table", 2))
            {
                ImGui::TableNextColumn();
                ImGui::Text("n");
                ImGui::TableNextColumn();
                ImGui::Text("error");
                for (int n = 0; n < max_points_number; ++n) {
                    ImGui::TableNextColumn();
                    ImGui::Text("%i", n+1);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.15f", lagrange_max_error[n]);
                }
                ImGui::TableNextColumn();
                ImGui::Text("------------------------------");
                ImGui::TableNextColumn();
                ImGui::TableNextColumn();
                ImGui::Text(("Newton error at " + std::to_string(lagrange_best_error_index)).c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.15f", newton_error);
                ImGui::TableNextColumn();
                ImGui::Text(("Chebyshev error at " + std::to_string(lagrange_best_error_index)).c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.15f", chebyshev_max_error);
                ImGui::TableNextColumn();
                ImGui::Text(("Spline error with n = " + std::to_string(spline_count)).c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.15f", spline_error);
                ImGui::TableNextColumn();
                ImGui::Text(("Best MNRP degree"));
                ImGui::TableNextColumn();
                ImGui::Text("%i", mnrp_best_degree);

                ImGui::EndTable();
            }

            ImGui::TableNextColumn();

            //if (ImPlot::BeginPlot("max_error")) {
            //    ImPlot::PlotBars("max_error", max_error.data(), 15);
            //    ImPlot::EndPlot();
            //}


            ImGui::EndTable();
        }

        ImGui::SliderInt("Degree of polynom", &current_plot_degree, 0, max_points_number - 1);
        if (ImGui::BeginTable("Interpolations", 2)) {
            ImGui::TableNextColumn();

            if (ImPlot::BeginPlot("Lagrange Interpolation")) {
                ImPlot::PlotLine("Ln(x)", test_xi, lagrange_test_yi[current_plot_degree], grid_size);
                ImPlot::PlotLine("f(x)", test_xi, target_yi, grid_size);
                ImPlot::PlotScatter("Grid points", lagrange_xi[current_plot_degree], lagrange_yi[current_plot_degree], current_plot_degree + 1);
                ImPlot::EndPlot();
            }

            ImGui::TableNextColumn();

            if (ImPlot::BeginPlot("Error of lagrange interpolation")) {
                ImPlot::PlotLine("Error at current degree", test_xi, lagrange_error[current_plot_degree], grid_size);
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Chebyshev roots")) {
                ImPlot::PlotLine("f(x)", test_xi, target_yi, grid_size);
                ImPlot::PlotLine("Lagrange on chebyshev roots", test_xi, chebyshev_test_yi, grid_size);
                ImPlot::PlotScatter("Grid points", chebyshev_grid_xi, chebyshev_grid_yi, lagrange_best_error_index + 1);
                ImPlot::EndPlot();
            }

            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Newton Interpolation")) {
                ImPlot::PlotLine("f(x)", test_xi, target_yi, grid_size);
                ImPlot::PlotLine("N_n0(x)", test_xi, newton_test_yi, grid_size);
                ImPlot::EndPlot();
            }

            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Trigonometry interpolation")) {
                ImPlot::PlotLine("f(x)", test_xi, target_yi, grid_size);
                ImPlot::PlotLine("trig(x)", test_xi, trigonometry_test_yi, grid_size);
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Trigonometry error")) {
                ImPlot::PlotLine("Error", test_xi, trigonometry_error_data, grid_size);
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Spline interpolation")) {
                ImPlot::PlotLine("f(x)", test_xi, target_yi, grid_size);
                ImPlot::PlotLine("S(x)", test_xi, spline_test_yi, grid_size);
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("Spline error")) {
                ImPlot::PlotLine("Error", test_xi, spline_error_data, grid_size);
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("MNRP")) {
                ImPlot::PlotLine("f(x)", test_xi, target_yi, grid_size);
                ImPlot::PlotLine("Pn", test_xi, P_yi, grid_size);
                ImPlot::PlotLine("Qm", test_xi, mnrp_result_yi, grid_size);
                ImPlot::EndPlot();
            }
            ImGui::TableNextColumn();
            if (ImPlot::BeginPlot("MNRP error")) {
                ImPlot::PlotLine("Error", test_xi, mnrp_error, grid_size);
                ImPlot::EndPlot();
            }
            ImGui::EndTable();

        }



       
        ImGui::End();

        window.clear();
        ImGui::SFML::Render(window);
        window.display();
    }


    ImGui::SFML::Shutdown();
    delete[] chebyshev_grid_xi;
    delete[] chebyshev_grid_yi;
    return 0;
}
