#include "kdv_5.h"
#include "service.h"
#include "gtest_ext.h"

#include <Eigen/Dense>

using Eigen::VectorXi;

TEST(ApproximationOrder, KawaharaEquation) {
    PRINTF("u_t + u u_x + u_xxx - u_xxxxx = 0\n");
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");

    const double a = 1.;
    const double b = 1.;
    const double c = 1.;
    const double d = 0.;

    KdV5Equation::coefficients eq_coeffs;

    eq_coeffs.c01 = a;
    eq_coeffs.c3 = b;
    eq_coeffs.c5 = -c;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;

    auto exact_sol = [a, b, c](double x, double t) -> double {
        const double arg = 1. / 2. * std::sqrt(b / 13. / c) * (x - 30. - 36. / 169. * b * b / c * t);
        return 105. / 169. * b * b / a / c / std::pow(std::cosh(arg), 4);
    };

    auto exact_sol_x = [exact_sol](double x, double t) -> double {
        const double arg = 1. / 2. / std::sqrt(13.) * (x - 10. - 36. / 169. * t);
        return -2. / std::sqrt(13) * std::tanh(arg) * exact_sol(x, t);
    };

    auto exact_sol_xx = [exact_sol](double x, double t) -> double {
        const double arg = 1. / 2. / std::sqrt(13.) * (x - 10. - 36. / 169. * t);
        return (5. / 13. * std::pow(std::tanh(arg), 2) - 1. / 13.) * exact_sol(x, t);
    };

    const std::vector<int> meshes{1000, 2000, 3000, 4000, 5000};
    const VectorXd steps_v = VectorXi::Map(meshes.data(), meshes.size()).cast<double>().cwiseInverse() * (x_max - x_min);
    std::vector<double> steps_stl_v(steps_v.size());
    VectorXd::Map(&steps_stl_v[0], steps_v.size()) = steps_v;

    std::ostringstream steps_str;
    steps_str << "steps: ";
    for (auto it = steps_stl_v.begin(); it < steps_stl_v.end() - 1; it++) {
        steps_str << *it << ", ";
    }
    steps_str << *(steps_stl_v.end() - 1) << "\n";

    PRINTF(steps_str.str().c_str());

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

    std::vector<double> errors;
    for (auto it = meshes.begin(); it < meshes.end(); it++) {
        KdV5Equation::meshProps mesh(x_min, x_max, t_max, *it);

        KdV5Equation kdv5(eq_coeffs, mesh,
                          iterations, 20, [exact_sol](double x) {
                              return exact_sol(x, 0.);
                          }, [exact_sol, x_min](double t) {
                              return exact_sol(x_min, t);
                          }, [exact_sol_x, x_min](double t) {
                              return exact_sol_x(x_min, t);
                          }, [exact_sol_xx, x_min](double t) {
                              return exact_sol_xx(x_min, t);
                          }, [exact_sol, x_max](double t) {
                              return exact_sol(x_max, t);
                          }, [exact_sol_x, x_max](double t) {
                              return exact_sol_x(x_max, t);
                          }, [exact_sol_xx, x_max](double t) {
                              return exact_sol_xx(x_max, t);
                          });
        kdv5.run();

        const ArrayXd mesh_x = kdv5.getXGrid();

        auto it2 = kdv5.getSolutions().rbegin();
        const double final_time = it2->first;
        const VectorXd final_solution = it2->second;

        VectorXd exact_solution(mesh_x.size());
        for (auto i = 0; i < mesh_x.size(); i++)
            exact_solution(i) = exact_sol(mesh_x(i), final_time);

        const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

        errors.push_back(error);
    }

    const VectorXd errors_v = VectorXd::Map(errors.data(), meshes.size());
    const ArrayXd steps_log = steps_v.array().log();
    const ArrayXd errors_log = errors_v.array().log();

    double lsf_a, lsf_b;
    std::tie(lsf_a, lsf_b) = computeLinearFitting(steps_log, errors_log);

    EXPECT_NEAR(lsf_a, 2., 8E-4);
    EXPECT_LT(std::exp(lsf_b), 1.);
}

TEST(ApproximationOrder, Kawahara6Equation) {
    PRINTF("u_t + u u_x + u_xxx - u_xxxxx = 0\n");
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");

    const double c01 = 2.;
    const double k = 0.25;
    const double c23 = -0.1;
    const double c12 = -10. * k * k * c23;
    const double c001 = 1.;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c01 = c01;
    eq_coeffs.c23 = c23;
    eq_coeffs.c12 = c12;
    eq_coeffs.c001 = c001;

    //std::cout << eq_coeffs.dump() << std::endl;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 60.;

    auto Q_func = [](double z) -> double {
        return 1./ (1. + std::exp(-z));
    };

    auto exact_sol = [c01, c001, c12, c23, k, Q_func](double x, double t) -> double {
        const double z = k * (x - 80.) - (144. * std::pow(k, 9) * c23 * c23 / c001 - 1. / 4. * k * c01 * c01 / c001) * t;
        return -12. * std::pow(k, 4) * c23 / c001 * (
                50. * std::pow(Q_func(z), 4)
                - 100. * std::pow(Q_func(z), 3)
                + 50. * std::pow(Q_func(z), 2)
                - 1.
        ) - c01 / 2. / c001;
    };

    auto exact_sol_x = [exact_sol, k](double x, double t) -> double {
        return 0.;
    };

    auto exact_sol_xx = [exact_sol, k](double x, double t) -> double {
        return 0.;
    };

    const std::vector<int> meshes{2000, 2500, 3000, 3500, 4000, 4500, 5000};
    const VectorXd steps_v = VectorXi::Map(meshes.data(), meshes.size()).cast<double>().cwiseInverse() * (x_max - x_min);
    std::vector<double> steps_stl_v(steps_v.size());
    VectorXd::Map(&steps_stl_v[0], steps_v.size()) = steps_v;

    std::ostringstream steps_str;
    steps_str << "steps: ";
    for (auto it = steps_stl_v.begin(); it < steps_stl_v.end() - 1; it++) {
        steps_str << *it << ", ";
    }
    steps_str << *(steps_stl_v.end() - 1) << "\n";

    PRINTF(steps_str.str().c_str());

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

    std::vector<double> errors;
    for (auto it = meshes.begin(); it < meshes.end(); it++) {
        KdV5Equation::meshProps mesh(x_min, x_max, t_max, *it);

        KdV5Equation kdv5(eq_coeffs, mesh,
                          iterations, 20, [exact_sol](double x) {
                              return exact_sol(x, 0.);
                          }, [exact_sol, x_min](double t) {
                              return exact_sol(x_min, t);
                          }, [exact_sol_x, x_min](double t) {
                              return exact_sol_x(x_min, t);
                          }, [exact_sol_xx, x_min](double t) {
                              return exact_sol_xx(x_min, t);
                          }, [exact_sol, x_max](double t) {
                              return exact_sol(x_max, t);
                          }, [exact_sol_x, x_max](double t) {
                              return exact_sol_x(x_max, t);
                          }, [exact_sol_xx, x_max](double t) {
                              return exact_sol_xx(x_max, t);
                          });
		//kdv5.setDebugStream(std::cout);
        kdv5.run();

        const ArrayXd mesh_x = kdv5.getXGrid();

        auto it2 = kdv5.getSolutions().rbegin();
        const double final_time = it2->first;
        const VectorXd final_solution = it2->second;

        VectorXd exact_solution(mesh_x.size());
        for (auto i = 0; i < mesh_x.size(); i++)
            exact_solution(i) = exact_sol(mesh_x(i), final_time);

        const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

        errors.push_back(error);
    }

    const VectorXd errors_v = VectorXd::Map(errors.data(), meshes.size());
    const ArrayXd steps_log = steps_v.array().log();
    const ArrayXd errors_log = errors_v.array().log();

    double lsf_a, lsf_b;
    std::tie(lsf_a, lsf_b) = computeLinearFitting(steps_log, errors_log);

    EXPECT_NEAR(lsf_a, 2., 3E-4);
    EXPECT_LT(std::exp(lsf_b), 1.);
}