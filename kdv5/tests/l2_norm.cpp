#include "kdv_5.h"
#include "gtest_ext.h"

std::tuple<double, double> compute_L2_norms(const std::map<double, VectorXd> &solution_list, double step) {
    VectorXd solution_initial = solution_list.begin()->second;
    VectorXd solution_finish = solution_list.rbegin()->second;
    double initial_norm = solution_initial.lpNorm<2>() * step;
    double finish_norm = solution_finish.lpNorm<2>() * step;

    return std::make_tuple(initial_norm, finish_norm);
}

TEST(L2norm, KawaharaEquation) {
    PRINTF("u_t + u u_x + u_xxx - u_xxxxx = 0\n"); // biswas
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");
    PRINTF("step: 0.02\n");

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
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

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

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-11);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 2E-5);
}

TEST(L2norm, KawaharaGeneral1Equation) {
    PRINTF("u_t + u^2 u_x + u_xxx - u_xxxxx = 0\n"); // 10.1016/j.physleta.2006.08.068
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");
    PRINTF("step: 0.02\n");

    const double a = 0.;
    const double b = 1.;
    const double c = 1.;
    const double d = 1.;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c3 = b;
    eq_coeffs.c5 = -c;
    eq_coeffs.c001 = d;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

    auto exact_sol = [b, c, d](double x, double t) -> double {
        const double arg = 1. / 2. * std::sqrt(1. / 5 * b / c) * (x - 30. - 4. / 25. * b * b / c * t);
        return -3. * b / std::sqrt(10. * d * c) / std::pow(std::cosh(arg), 2);
    };

    auto exact_sol_x = [exact_sol, b, c](double x, double t) -> double {
        const double arg = 1. / 2. * std::sqrt(1. / 5 * b / c) * (x - 30. - 4. / 25. * b * b / c * t);
        return -std::sqrt(1. / 5. * b / c) * std::tanh(arg) * exact_sol(x, t);
    };

    auto exact_sol_xx = [exact_sol, b, c](double x, double t) -> double {
        const double arg = 1. / 2. * std::sqrt(1. / 5 * b / c) * (x - 30. - 4. / 25. * b * b / c * t);
        return 1. / 10. * b / c * (3. * std::pow(std::tanh(arg), 2) - 1.) * exact_sol(x, t);
    };

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-10);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 2E-5);
}

TEST(L2norm, KawaharaGeneral2Equation) {
    PRINTF("u_t + u u_x + u^2 u_x + u_xxx - u_xxxxx = 0\n"); // 10.1016/j.physleta.2006.08.068
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");
    PRINTF("step: 0.02\n");

    const double a = 1.;
    const double b = 1.;
    const double c = 1.;
    const double d = 1.;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c01 = a;
    eq_coeffs.c3 = b;
    eq_coeffs.c5 = -c;
    eq_coeffs.c001 = d;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

    auto exact_sol = [a, b, c, d](double x, double t) -> double {
        const double arg = std::sqrt(5.) / 10. * std::sqrt(b / c) * (x - 60. + 1. / 4. * a * a / d * t - 4. / 25. * b * b / c * t);
        return -3. / std::sqrt(10.) * b / std::sqrt(c * d) / std::pow(std::cosh(arg), 2) - 1. / 2. * a / d;
    };

    auto exact_sol_x = [a, b, c, d](double x, double t) -> double {
        const double arg = std::sqrt(5.) / 10. * std::sqrt(b / c) * (x - 60. + 1. / 4. * a * a / d * t - 4. / 25. * b * b / c * t);
        return 3. / 10. * std::sqrt(2.) * b / c * std::sqrt(b / d) * std::sinh(arg) / std::pow(std::cosh(arg), 3);
    };

    auto exact_sol_xx = [a, b, c, d](double x, double t) -> double {
        const double arg = std::sqrt(5.) / 10. * std::sqrt(b / c) * (x - 60. + 1. / 4. * a * a / d * t - 4. / 25. * b * b / c * t);
        return -3. / 100. * std::sqrt(10.) * b * b / c / std::sqrt(c * d) * (2. * std::pow(std::cosh(arg), 2) - 3.) / std::pow(std::cosh(arg), 4);
    };

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-7);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 2E-5);
}

TEST(L2norm, KawaharaGeneral3Equation) {
    PRINTF("u_t + u_x u_xx + u_xxxxx = 0\n"); // 10.4208/jpde.v25.n4.4
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");
    PRINTF("step: 0.02\n");

    const double a = 0.;
    const double b = 0.;
    const double c = 1.;
    const double c12 = 1.;
    const double d = 0.;
    const double k = 0.25;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c01 = a;
    eq_coeffs.c3 = b;
    eq_coeffs.c5 = c;
    eq_coeffs.c001 = d;
    eq_coeffs.c12 = c12;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

    auto exact_sol = [d, c12, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return 60. * std::pow(k, 2) / c12 * (1. - std::pow(std::tanh(arg), 2));
    };

    auto exact_sol_x = [exact_sol, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return -2. * k * std::tanh(arg) * exact_sol(x, t);
    };

    auto exact_sol_xx = [exact_sol, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return 2. * std::pow(k, 2) * (3. * std::pow(std::tanh(arg), 2) - 1) * exact_sol(x, t);
    };

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-4);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 8E-5);
}

TEST(L2norm, KawaharaGeneral4Equation) {
    PRINTF("u_t + u^2 u_x + u * u_xxx + u_xxxxx = 0\n"); // 10.4208/jpde.v25.n4.4
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");
    PRINTF("step: 0.02\n");

    const double a = 0.;
    const double b = 0.;
    const double c = 1.;
    const double c12 = 0.;
    const double c03 = 1.;
    const double d = 1. / 10. * c03 * c03;
    const double k = 0.25;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c01 = a;
    eq_coeffs.c3 = b;
    eq_coeffs.c5 = c;
    eq_coeffs.c001 = d;
    eq_coeffs.c12 = c12;
    eq_coeffs.c03 = c03;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

    auto exact_sol = [d, c12, c03, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return 60. * std::pow(k, 2) / c03 * (1. - std::pow(std::tanh(arg), 2));
    };

    auto exact_sol_x = [exact_sol, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return -2. * k * std::tanh(arg) * exact_sol(x, t);
    };

    auto exact_sol_xx = [exact_sol, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return 2. * std::pow(k, 2) * (3. * std::pow(std::tanh(arg), 2) - 1) * exact_sol(x, t);
    };

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-4);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 3E-4);
}

TEST(L2norm, KawaharaGeneral5Equation) {
    PRINTF("u_t + u^2 u_x + u * u_xxx + u_xxxxx = 0\n"); // 10.4208/jpde.v25.n4.4
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:100]\n");
    PRINTF("step: 0.02\n");

    const double c = 1.;
    const double c12 = 0.5;
    const double c03 = 0.5;
    const double d = 1. / 10. * (c12 + c03) * c03;
    const double k = 0.25;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c03 = c03;
    eq_coeffs.c12 = c12;
    eq_coeffs.c5 = c;
    eq_coeffs.c001 = d;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

    auto exact_sol = [d, c12, c03, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return 60. * std::pow(k, 2) / (c12 + c03) * (1. - std::pow(std::tanh(arg), 2));
    };

    auto exact_sol_x = [exact_sol, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return -2. * k * std::tanh(arg) * exact_sol(x, t);
    };

    auto exact_sol_xx = [exact_sol, k](double x, double t) -> double {
        const double arg = k * (x - 30.) - 16. * std::pow(k, 5) * t;
        return 2. * std::pow(k, 2) * (3. * std::pow(std::tanh(arg), 2) - 1) * exact_sol(x, t);
    };

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-5);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 3E-4);
}

TEST(L2norm, KawaharaGeneral6Equation) {
    PRINTF("u_t + 2 * u u_x + u^2 u_x + 0.0625 * u_x u_xx - 0.1 * u_xx u_xxx = 0\n"); // maple
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:50]\n");
    PRINTF("step: 0.02\n");

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

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 50.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size);

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

    KdV5Equation::iterationsProps iterations = {50, 1E-16};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-6);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 2E-6);
}

TEST(L2norm, KawaharaGeneral7Equation) {
    PRINTF("u_t + 0.158567 * u^2 u_x + 54.4794 * u_x u_xx - 106.196 * u_xx u_xxx + 60.7178 * u_x u_xxxx + u u_xxxxx = 0\n"); // maple
    PRINTF("x: [0:100]\n");
    PRINTF("t: [0:10]\n");
    PRINTF("step: 0.02\n");

    const double k = 0.23;
    const double omega = -0.10844852666666666667e-2;
    const double a0 = 0;
    const double a1 = 1;
    const double a2 = -1;

    KdV5Equation::coefficients eq_coeffs;
    eq_coeffs.c5 = 1.;
    eq_coeffs.c03 = 1.;
    eq_coeffs.c12 = 1.;

    eq_coeffs.c001 = -0.5522760e-1;
    eq_coeffs.c3 = -.14203333333333333333;
    eq_coeffs.c01 = -0.28058160000000000000e-1;

    const double x_min = 0.;
    const double x_max = 100.;
    const double t_max = 100.;
    unsigned int mesh_size = 5000;
    KdV5Equation::meshProps mesh(x_min, x_max, t_max, mesh_size, 1E0);

    auto Q_func = [](double z) -> double {
        return 1. / (1. + std::exp(-z));
    };

    auto exact_sol = [a0, a1, a2, k, omega, Q_func](double x, double t) -> double {
        const double z = k * (x - 60.) - omega * t;
        return a0 + a1 * std::pow(Q_func(z), 1) + a2 * std::pow(Q_func(z), 2);
    };

    auto exact_sol_x = [](double x, double t) -> double {
        return 0.;
    };

    auto exact_sol_xx = [](double x, double t) -> double {
        return 0.;
    };

    KdV5Equation::iterationsProps iterations = {100, 1E-16, 100};

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
    //std::cout << kdv5.dump() << std::endl;
    kdv5.run();

    double initial_norm, finish_norm;
    std::tie(initial_norm, finish_norm) = compute_L2_norms(kdv5.getSolutions(), mesh.x_step);

    PRINTF("initial norm: ");
    //std::cout << std::fabs(initial_norm) << "\n";

    PRINTF("|initial norm - finish norm|: ");
    //std::cout << std::fabs(finish_norm - initial_norm) << "\n";
    PRINTF("|initial norm - finish norm| / initial norm: ");
    //std::cout << std::fabs(finish_norm - initial_norm) / initial_norm << "\n";
    EXPECT_NEAR(initial_norm, finish_norm, 1E-6);


    const ArrayXd mesh_x = kdv5.getXGrid();

    auto it = kdv5.getSolutions().rbegin();
    const double final_time = it->first;
    const VectorXd final_solution = it->second;

    VectorXd exact_solution(mesh_x.size());
    for (auto i = 0; i < mesh_x.size(); i++)
        exact_solution(i) = exact_sol(mesh_x(i), final_time);

    const double error = (final_solution - exact_solution).cwiseAbs().maxCoeff();

    EXPECT_NEAR(error, 0, 9E-5);
}