#include "kdv_5.h"

#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

int main() {
    try {
        KdV5Equation::modelParameters model_params;
        model_params.a = 1.;
        model_params.l = 800.;
        model_params.h_0 = 50.;
        model_params.h_i = 3;
        const KdV5Equation::expParameters small_params = model_params;
        const KdV5Equation::coefficients coeffs = small_params;

        std::cout
                << "[model parameters]\n" << model_params << "\n"
                << "[small parameters]\n" << small_params << "\n"
                << "[equation coefficients]\n" << coeffs << "\n";

        const double x_min = 0.;
        const double x_max = 100.;
        const double t_max = 100.;

        auto u_initial = [](double x) -> double {
            return 1. / std::pow(std::cosh(x - 50.), 2);
        };

        auto zero_func = [](double x) -> double {
            return 0.;
        };

        const auto total_start = Clock::now();
		KdV5Equation::meshProps mesh(x_min, x_max, t_max, 5000, 1E0);
		KdV5Equation::iterationsProps iterations = {100, 1E-16, 100};

		const int plot_count = 10;
		KdV5Equation kdv5(coeffs, mesh, iterations, plot_count,
						  u_initial, zero_func, zero_func, zero_func, zero_func, zero_func, zero_func);
//      kdv5.setDebugStream(std::cout);
		kdv5.setSolutionFile("solution_example_1.data");
		std::cout << kdv5 << std::endl;

		std::cout << "Calculation start..." << std::endl;
		auto t1 = Clock::now();
		kdv5.run();
		auto t2 = Clock::now();
		std::cout << "Calculation done" << std::endl;

		std::cout << "Total time: "
				  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Exception:" << std::endl << ex.what() << std::endl;
    }

	std::cout << "Press any key..." << std::endl;
    std::getchar();

    return 0;
}
