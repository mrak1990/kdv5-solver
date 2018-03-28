#include "service.h"
#include "kdv_5.h"

#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

int main() {
    try {
        KdV5Equation::modelParameters model_params_i;
        model_params_i.a = 1.;
        model_params_i.l = 800.;
        model_params_i.h_0 = 50.;
        model_params_i.h_i = 3;
        const KdV5Equation::expParameters small_params_i = model_params_i;
        const KdV5Equation::coefficients coeffs_i = small_params_i;

        KdV5Equation::modelParameters model_params_w = model_params_i;
        model_params_w.h_i = 0;
        const KdV5Equation::expParameters small_params_w = model_params_w;
        const KdV5Equation::coefficients coeffs_w = small_params_w;

        std::cout
                << "[model parameters]\n" << model_params_i << "\n"
                << "[small parameters (ice)]\n" << small_params_i << "\n"
                << "[small parameters (water)]\n" << small_params_w << "\n"
                << "[equation coefficients (ice)]\n" << coeffs_i << "\n"
                << "[equation coefficients (water)]\n" << coeffs_w << std::endl;

        const double x_min = 0.;
        const double x_max = 100.;
        const double t_max = 100.;

        auto u_initial = [](double x) -> double {
            return 1. / std::pow(std::cosh(x - 50.), 2);
        };

        auto zero_func = [](double x) -> double {
            return 0.;
        };

        std::vector<int> meshes{5000};
        std::vector<double> errors;

        const auto total_start = Clock::now();
        for (auto it = meshes.begin(); it < meshes.end(); it++) {
            KdV5Equation::meshProps mesh(x_min, x_max, t_max, *it, 1E0);
            KdV5Equation::iterationsProps iterations = {100, 1E-16, 100};

            const int plot_count = 10;
            KdV5Equation kdv5(coeffs_i, coeffs_w, mesh, iterations, plot_count,
                              u_initial, zero_func, zero_func, zero_func, zero_func, zero_func, zero_func);
//            kdv5.setDebugStream(std::cout);
            kdv5.setSolutionFile("solution.data");
            std::cout << kdv5 << std::endl;

            std::cout << "Calculation start..." << std::endl;
            auto t1 = Clock::now();
            kdv5.run();
            auto t2 = Clock::now();
            std::cout << "Calculation done" << std::endl;

            std::cout << "Total time: "
                      << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Exception:" << std::endl << ex.what() << std::endl;
    }

    std::getchar();

    return 0;
}
