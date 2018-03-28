#include <fstream>

#include "kdv_5.h"

using Eigen::Sequential;

KdV5Equation::KdV5Equation(coefficients coefficients,
                           meshProps mesh2, iterationsProps iterations, unsigned int plot_count,
                           func initial, func u_a, func ux_a,
                           func uxx_a, func u_b, func ux_b, func uxx_b) :
        coeffs(coefficients), mesh(mesh2), iterations(iterations),
        plot_count(plot_count),
        initial(initial), u_a(u_a), ux_a(ux_a), uxx_a(uxx_a), u_b(u_b), ux_b(ux_b), uxx_b(uxx_b),
        solver(mesh2.x_size + 1) {
    solution_u.resize(mesh2.x_size + 1);
    old_solution_u.resize(mesh2.x_size + 1);
    system_matrix.resize(10, mesh2.x_size + 1);
    system_rhs.resize(mesh2.x_size + 1);

    out.timestep_per_plot = static_cast<unsigned int>((mesh2.t_max + mesh2.t_step) / mesh2.t_step / plot_count);
    out.timestep_number_to_plot = 0;

    _x_step_p_m1 = 1. / mesh2.x_step;
    _x_step_p_m2 = 1. / std::pow(mesh2.x_step, 2);
    _x_step_p_m3 = 1. / std::pow(mesh2.x_step, 3);
    _x_step_p_m4 = 1. / std::pow(mesh2.x_step, 4);
    _x_step_p_m5 = 1. / std::pow(mesh2.x_step, 5);
    _t_step_p_m1 = 1. / mesh2.t_step;
}

void KdV5Equation::run() {
    if (!out.solution_file.empty()) {
        std::ofstream ofs;
        ofs.open(out.solution_file, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
    }

    xGrid = VectorXd::LinSpaced(Sequential, mesh.x_size + 1, mesh.x_min, mesh.x_max);

    for (auto i = 0; i <= mesh.x_size; i++)
        old_solution_u(i) = initial(xGrid(i));
    solution_u = old_solution_u;

    double solution_norm_1 = solution_u.lpNorm<1>();
    double solution_norm_2 = solution_u.lpNorm<2>();

    if (out.debug_stream)
        *out.debug_stream
                << "L1 initial norm: " << mesh.x_step * solution_norm_1 << std::endl
                << "l1 initial norm: " << solution_norm_1 << std::endl
                << "L2 initial norm: " << mesh.x_step * solution_norm_2 << std::endl
                << "l2 initial norm: " << solution_norm_2 << std::endl << std::endl;

    time = 0;
    unsigned int timestep_counter = 0;
    dumpSolution();
    while (timestep_counter++, time = mesh.t_step * timestep_counter, time <= mesh.t_max + mesh.t_step) {
        if (out.debug_stream)
            *out.debug_stream << "t=" << time << "\n";

        bool convergence = false;
        unsigned int cycle = 1;
        while (cycle++, cycle <= iterations.max) {
            assemble_system();
            assemble_rhs();
            solver.solve(system_matrix.data(), system_rhs.data());
            solution_u -= system_rhs;

            double error = system_rhs.lpNorm<2>();
            if (error / solution_norm_1 < iterations.max_error) {
                convergence = true;
                break;
            }

            if (iterations.exact > 0 && cycle == iterations.exact) {
                convergence = true;
                break;
            }
        }

        if (!convergence)
            throw std::logic_error("Iterations exceeds");

        old_solution_u = solution_u;

        if (timestep_counter == out.timestep_number_to_plot) {
            dumpSolution();
        }
    }

    solution_norm_1 = solution_u.lpNorm<1>();
    solution_norm_2 = solution_u.lpNorm<2>();

    if (out.debug_stream)
        *out.debug_stream
                << "L1 final norm: " << mesh.x_step * solution_norm_1 << std::endl
                << "l1 final norm: " << solution_norm_1 << std::endl
                << "L2 final norm: " << mesh.x_step * solution_norm_2 << std::endl
                << "l2 final norm: " << solution_norm_2 << std::endl << std::endl;
}

void KdV5Equation::setDebugStream(std::ostream &stream) {
    out.debug_stream = &stream;
}

void KdV5Equation::setInstantDump(bool instant) {
    out.instant_dump = instant;
}

void KdV5Equation::setSolutionFile(std::string file) {
    out.solution_file = file;
}

const ArrayXd &KdV5Equation::getXGrid() const {
    return xGrid;
}

const std::map<double, VectorXd> &KdV5Equation::getSolutions() const {
    return solution_list;
}

void KdV5Equation::dumpSolution() {
    if (out.instant_dump && !out.solution_file.empty()) {
        std::ofstream ofs(out.solution_file, std::ofstream::out | std::ofstream::app);
		std::cout << "Write to file (t=" << time << ")... ";

        MatrixXd m(xGrid.size(), 3);
        m.col(0).setConstant(time);
        m.col(1) = xGrid;
        m.col(2) = solution_u;

        ofs << m << std::endl << std::endl << std::endl;
        ofs.close();
        std::cout << "done" << std::endl;
    }

    solution_list.insert(std::pair<double, VectorXd>(time, old_solution_u));
    out.timestep_number_to_plot += out.timestep_per_plot;
}

void KdV5Equation::assemble_system() {
    VectorXd &u = solution_u;
    VectorXd &u_p = old_solution_u;
    const auto &n = mesh.x_size;

    const coefficients &c_i = coeffs;
    const double val1 = c_i.c5 / 4. * _x_step_p_m5;
    const double val2 = c_i.c3 / 4. * _x_step_p_m3 - c_i.c5 * _x_step_p_m5;

    // J_j-3_j
    system_matrix.row(3).setConstant(val1);
    // J_0_3
    system_matrix(3, 3) = 0.;
    // J_1_4
    system_matrix(3, 4) = 0.;
    // J_2_5
    system_matrix(3, 5) = 0.;

    // J_j-2_j
    system_matrix.row(4).setConstant(val2);
    // J_0_2
    system_matrix(4, 2) = 0.;
    // J_1_3
    system_matrix(4, 3) = 0.;
    // J_2_4
    system_matrix(4, 4) = 0.;
    // J_n-2_n
    system_matrix(4, n) = 2. * _x_step_p_m2;

    // J_0_1
    system_matrix(5, 1) = 0.;
    // J_1_2
    system_matrix(5, 2) = -1. / 2. * _x_step_p_m1;
    // J_2_3
    system_matrix(5, 3) = -_x_step_p_m2;
    // J_j-1_j
    for (unsigned int j = 4; j <= n - 2; j++)
        system_matrix(5, j) =
                c_i.c01 / 8. * _x_step_p_m1
                * (
                        u(j - 1)
                        + u_p(j - 1)
                )
                - c_i.c3 / 2. * _x_step_p_m3
                + c_i.c5 * 5. / 4. * _x_step_p_m5
                + c_i.c001 / 16. * _x_step_p_m1
                  * (
                          u(j - 1)
                          + u_p(j - 1)
                  )
                  * (
                          u(j - 1)
                          + u_p(j - 1)
                  )
                - c_i.c03 / 4. * _x_step_p_m3
                  * (
                          u(j - 1)
                          + u_p(j - 1)
                  );
    // J_n-2_n-1
    system_matrix(5, n - 1) = -5. * _x_step_p_m2;
    // J_n-1_n
    system_matrix(5, n) = 3. / 2. * _x_step_p_m1;

    // J_0_0
    system_matrix(6, 0) = 1.;
    // J_1_1
    system_matrix(6, 1) = 2. * _x_step_p_m1;
    // J_2_2
    system_matrix(6, 2) = 4. * _x_step_p_m2;
    // J_j_j
    for (unsigned int j = 3; j <= n - 3; j++)
        system_matrix(6, j) =
                _t_step_p_m1
                + c_i.c01 / 8. * _x_step_p_m1
                  * (
                          u(j + 1) - u(j - 1)
                          + u_p(j + 1) - u_p(j - 1)
                  )
                + c_i.c001 / 8. * _x_step_p_m1
                  * (
                          u(j)
                          + u_p(j)
                  )
                  * (
                          u(j + 1) - u(j - 1)
                          + u_p(j + 1) - u_p(j - 1)
                  )
                - c_i.c12 / 4. * _x_step_p_m3
                  * (
                          u(j + 1) - u(j - 1)
                          + u_p(j + 1) - u_p(j - 1)
                  )
                + c_i.c03 / 8. * _x_step_p_m3
                  * (
                          -u(j - 2) + 2. * u(j - 1) - 2. * u(j + 1) + u(j + 2)
                          - u_p(j - 2) + 2. * u_p(j - 1) - 2. * u_p(j + 1) + u_p(j + 2)
                  )
                - c_i.c23 / 4. * _x_step_p_m5
                  * (
                          -u(j - 2) + 2. * u(j - 1) - 2. * u(j + 1) + u(j + 2)
                          - u_p(j - 2) + 2. * u_p(j - 1) - 2. * u_p(j + 1) + u_p(j + 2)
                  )
                + c_i.c14 * 3. / 4. * _x_step_p_m5
                  * (
                          u(j + 1) - u(j - 1)
                          + u_p(j + 1) - u_p(j - 1)
                  )
                + c_i.c4 * 3. * _x_step_p_m4
                - c_i.c2 * _x_step_p_m2
                + c_i.c05 / 8. * _x_step_p_m5
                  * (
                          -u(j - 3) + 4. * u(j - 2) - 5. * u(j - 1) + 5. * u(j + 1) - 4. * u(j + 2) + u(j + 3)
                          - u_p(j - 3) + 4. * u_p(j - 2) - 5. * u_p(j - 1) + 5. * u_p(j + 1) - 4. * u_p(j + 2) +
                          u_p(j + 3)
                  )
                - c_i.c22 * _x_step_p_m4
                  * (
                          u(j - 1) - 2. * u(j) + u(j + 1)
                          + u_p(j - 1) - 2. * u_p(j) + u_p(j + 1)
                  )
                + c_i.c04 / 4. * _x_step_p_m4
                  * (
                          u(j - 2) - 4. * u(j - 1) + 6. * u(j) - 4. * u(j + 1) + u(j + 2)
                          + u_p(j - 2) - 4. * u_p(j - 1) + 6. * u_p(j) - 4. * u_p(j + 1) + u_p(j + 2)
                          + 6.
                            * (
                                    u(j)
                                    + u_p(j)
                            )
                  );
    // J_n-2_n-2
    system_matrix(6, n - 2) = 4. * _x_step_p_m2;
    // J_n-1_n-1
    system_matrix(6, n - 1) = -2. * _x_step_p_m1;
    // J_n_n
    system_matrix(6, n) = 1.;

    // J_1_0
    system_matrix(7, 0) = -3. / 2 * _x_step_p_m1;
    // J_2_1
    system_matrix(7, 1) = -5. * _x_step_p_m2;
    // J_j+1_j
    for (unsigned int j = 2; j <= n - 4; j++)
        system_matrix(7, j) = -system_matrix(5, j + 2);
    // J_n-2_n-3
    system_matrix(7, n - 3) = -_x_step_p_m2;
    // J_n-1_n-2
    system_matrix(7, n - 2) = 1. / 2. * _x_step_p_m1;
    // J_n_n-1
    system_matrix(7, n - 1) = 0.;

    // J_j+2_j
    system_matrix.row(8).setConstant(-val2);
    // J_2_0
    system_matrix(8, 0) = 2. * _x_step_p_m2;
    // J_n-2_n-4
    system_matrix(8, n - 4) = 0.;
    // J_n-1_n-3
    system_matrix(8, n - 3) = 0.;
    // J_n_n-2
    system_matrix(8, n - 2) = 0.;

    // J_j+3_j
    system_matrix.row(9).setConstant(-val1);
    // J_n-2_n-5
    system_matrix(9, n - 5) = 0.;
    // J_n-1_n-4
    system_matrix(9, n - 4) = 0.;
    // J_n_n-3
    system_matrix(9, n - 3) = 0.;

    // J_j-1_j
    for (unsigned int j = 4; j <= n - 2; j++)
        system_matrix(5, j) +=
                c_i.c12 / 8. * _x_step_p_m3
                * (
                        u(j - 2) - 2. * u(j - 1) + u(j)
                        + u_p(j - 2) - 2. * u_p(j - 1) + u_p(j)
                        + u(j) - u(j - 2)
                        + u_p(j) - u_p(j - 2)
                )
                + c_i.c23 / 8. * _x_step_p_m5
                  * (
                          -u(j - 3) + 2. * u(j - 2) - 2. * u(j) + u(j + 1)
                          - u_p(j - 3) + 2. * u_p(j - 2) - 2. * u_p(j) + u_p(j + 1)
                          - 2.
                            * (
                                    u(j - 2) - 2. * u(j - 1) + u(j)
                                    + u_p(j - 2) - 2. * u_p(j - 1) + u_p(j)
                            )
                  )
                + c_i.c14 / 8. * _x_step_p_m5
                  * (
                          u(j - 3) - 4. * u(j - 2) + 6. * u(j - 1) - 4. * u(j) + u(j + 1)
                          + u_p(j - 3) - 4. * u_p(j - 2) + 6. * u_p(j - 1) - 4. * u_p(j) + u_p(j + 1)
                          - 4. * (
                                  u(j) - u(j - 2)
                                  + u_p(j) - u_p(j - 2)
                          )
                  )
                - c_i.c4 * 2. * _x_step_p_m4
                + c_i.c2 / 2. * _x_step_p_m2
                + c_i.c05 * 5. / 8. * _x_step_p_m5
                  * (
                          u(j - 1)
                          + u_p(j - 1)
                  )
                + c_i.c13 / 16. * _x_step_p_m4
                  * (
                          -u(j - 3) + 2. * u(j - 2) - 2. * u(j) + u(j + 1)
                          - u_p(j - 3) + 2. * u_p(j - 2) - 2. * u_p(j) + u_p(j + 1)
                          - 2.
                            * (
                                    u(j) - u(j - 2)
                                    + u_p(j) - u_p(j - 2)
                            )
                  )
                + c_i.c22 / 2. * _x_step_p_m4
                  * (
                          u(j - 2) - 2. * u(j - 1) + u(j)
                          + u_p(j - 2) - 2. * u_p(j - 1) + u_p(j)
                  )
                - c_i.c04 * _x_step_p_m4
                  * (
                          u(j - 1)
                          + u_p(j - 1)
                  );

    // J_j+1_j
    for (unsigned int j = 2; j <= n - 4; j++)
        system_matrix(7, j) +=
                c_i.c12 / 8. * _x_step_p_m3
                * (
                        -1. * (
                                u(j) - 2. * u(j + 1) + u(j + 2)
                                + u_p(j) - 2. * u_p(j + 1) + u_p(j + 2)

                        )
                        + u(j + 2) - u(j)
                        + u_p(j + 2) - u_p(j)
                )
                + c_i.c23 / 8. * _x_step_p_m5
                  * (
                          -u(j - 1) + 2. * u(j) - 2. * u(j + 2) + u(j + 3)
                          - u_p(j - 1) + 2. * u_p(j) - 2. * u_p(j + 2) + u_p(j + 3)
                          + 2.
                            * (
                                    u(j) - 2. * u(j + 1) + u(j + 2)
                                    + u_p(j) - 2. * u_p(j + 1) + u_p(j + 2)
                            )
                  )
                - c_i.c14 / 8. * _x_step_p_m5
                  * (
                          u(j - 1) - 4. * u(j) + 6. * u(j + 1) - 4. * u(j + 2) + u(j + 3)
                          + u_p(j - 1) - 4. * u_p(j) + 6. * u_p(j + 1) - 4. * u_p(j + 2) + u_p(j + 3)
                          + 4. * (
                                  u(j + 2) - u(j)
                                  + u_p(j + 2) - u_p(j)
                          )
                  )
                - c_i.c4 * 2. * _x_step_p_m4
                + c_i.c2 / 2. * _x_step_p_m2
                - c_i.c05 * 5. / 8. * _x_step_p_m5
                  * (
                          u(j + 1)
                          + u_p(j + 1)
                  )
                + c_i.c13 / 16. * _x_step_p_m4
                  * (
                          -1.
                          * (
                                  -u(j - 1) + 2. * u(j) - 2. * u(j + 2) + u(j + 3)
                                  - u_p(j - 1) + 2. * u_p(j) - 2. * u_p(j + 2) + u_p(j + 3)
                          )
                          + 2.
                            * (
                                    u(j + 2) - u(j)
                                    + u_p(j + 2) - u_p(j)
                            )
                  )
                + c_i.c22 / 2. * _x_step_p_m4
                  * (
                          u(j) - 2. * u(j + 1) + u(j + 2)
                          + u_p(j) - 2. * u_p(j + 1) + u_p(j + 2)
                  )
                - c_i.c04 * _x_step_p_m4
                  * (
                          u(j + 1)
                          + u_p(j + 1)
                  );

    // J_j-3_j
    for (unsigned int j = 6; j <= n; j++)
        system_matrix(3, j) +=
                +c_i.c05 / 8. * _x_step_p_m5
                * (
                        u(j - 3)
                        + u_p(j - 3)
                );
    // J_j+3_j
    for (unsigned int j = 0; j <= n - 5; j++)
        system_matrix(9, j) +=
                -c_i.c05 / 8. * _x_step_p_m5
                * (
                        u(j + 3)
                        + u_p(j + 3)
                );


    // -----------------------
    // J_j-2_j
    for (unsigned int j = 5; j <= n - 1; j++)
        system_matrix(4, j) +=
                c_i.c03 / 8. * _x_step_p_m3
                * (
                        u(j - 2)
                        + u_p(j - 2)
                )
                - c_i.c05 / 2. * _x_step_p_m5
                  * (
                          u(j - 2)
                          + u_p(j - 2)
                  )
                + c_i.c14 / 8. * _x_step_p_m5
                  * (
                          u(j - 1) - u(j - 3)
                          + u_p(j - 1) - u_p(j - 3)
                  )
                + c_i.c23 / 8. * _x_step_p_m5
                  * (
                          u(j - 3) - 2. * u(j - 2) + u(j - 1)
                          + u_p(j - 3) - 2. * u_p(j - 2) + u_p(j - 1)
                  )
                + c_i.c4 / 2. * _x_step_p_m4
                + c_i.c04 / 4. * _x_step_p_m4
                  * (
                          u(j - 2)
                          + u_p(j - 2)
                  )
                + c_i.c13 / 16. * _x_step_p_m4
                  * (
                          u(j - 1) - u(j - 3)
                          + u_p(j - 1) - u_p(j - 3)
                  );
    // J_j+2_j
    for (unsigned int j = 1; j <= n - 5; j++)
        system_matrix(8, j) +=
                -c_i.c03 / 8. * _x_step_p_m3
                * (
                        u(j + 2)
                        + u_p(j + 2)
                )
                + c_i.c05 / 2. * _x_step_p_m5
                  * (
                          u(j + 2)
                          + u_p(j + 2)
                  )
                + c_i.c14 / 8. * _x_step_p_m5
                  * (
                          u(j + 3) - u(j + 1)
                          + u_p(j + 3) - u_p(j + 1)
                  )
                - c_i.c23 / 8. * _x_step_p_m5
                  * (
                          u(j + 1) - 2. * u(j + 2) + u(j + 3)
                          + u_p(j + 1) - 2. * u_p(j + 2) + u_p(j + 3)
                  )
                + c_i.c4 / 2. * _x_step_p_m4
                + c_i.c04 / 4. * _x_step_p_m4
                  * (
                          u(j + 2)
                          + u_p(j + 2)
                  )
                - c_i.c13 / 16. * _x_step_p_m4
                  * (
                          u(j + 3) - u(j + 1)
                          + u_p(j + 3) - u_p(j + 1)
                  );
}

void KdV5Equation::assemble_rhs() {
    VectorXd &u = solution_u;
    VectorXd &u_p = old_solution_u;
    const auto &n = mesh.x_size;
    const coefficients &c = coeffs;

    system_rhs(0) = u(0) - u_a(time);
    system_rhs(1) = (-3. * u(0) + 4. * u(1) - u(2)) / 2. * _x_step_p_m1 - ux_a(time);
    system_rhs(2) = (2. * u(0) - 5. * u(1) + 4. * u(2) - u(3)) * _x_step_p_m2 - uxx_a(time);
    for (unsigned int i = 3; i <= n - 3; i++) {
        system_rhs(i) =
                    (u(i) - u_p(i)) * _t_step_p_m1
                    + c.c3 / 4. * _x_step_p_m3
                      * (
                              -u(i - 2) + 2. * u(i - 1) - 2. * u(i + 1) + u(i + 2)
                              - u_p(i - 2) + 2. * u_p(i - 1) - 2. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c5 / 4. * _x_step_p_m5
                      * (
                              -u(i - 3) + 4. * u(i - 2) - 5. * u(i - 1) + 5. * u(i + 1) - 4. * u(i + 2) + u(i + 3)
                              - u_p(i - 3) + 4. * u_p(i - 2) - 5. * u_p(i - 1) + 5. * u_p(i + 1) - 4. * u_p(i + 2) +
                              u_p(i + 3)
                      )
                    + c.c01 / 8. * _x_step_p_m1
                      * (
                              u(i)
                              + u_p(i)
                      )
                      * (
                              u(i + 1) - u(i - 1)
                              + u_p(i + 1) - u_p(i - 1)
                      )
                    + c.c03 / 8. * _x_step_p_m3
                      * (u(i) + u_p(i))
                      * (
                              -u(i - 2) + 2. * u(i - 1) - 2. * u(i + 1) + u(i + 2)
                              - u_p(i - 2) + 2. * u_p(i - 1) - 2. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c12 / 8. * _x_step_p_m3
                      * (
                              u(i + 1) - u(i - 1)
                              + u_p(i + 1) - u_p(i - 1)
                      )
                      * (
                              u(i - 1) - 2. * u(i) + u(i + 1)
                              + u_p(i - 1) - 2. * u_p(i) + u_p(i + 1)
                      )
                    + c.c05 / 8. * _x_step_p_m5
                      * (
                              u(i)
                              + u_p(i)
                      )
                      * (
                              -u(i - 3) + 4. * u(i - 2) - 5. * u(i - 1) + 5. * u(i + 1) - 4. * u(i + 2) + u(i + 3)
                              - u_p(i - 3) + 4. * u_p(i - 2) - 5. * u_p(i - 1) + 5. * u_p(i + 1) - 4. * u_p(i + 2) +
                              u_p(i + 3)
                      )
                    + c.c14 / 8. * _x_step_p_m5
                      * (
                              u(i + 1) - u(i - 1)
                              + u_p(i + 1) - u_p(i - 1)
                      )
                      * (
                              u(i - 2) - 4. * u(i - 1) + 6. * u(i) - 4. * u(i + 1) + u(i + 2)
                              + u_p(i - 2) - 4. * u_p(i - 1) + 6. * u_p(i) - 4. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c23 / 8. * _x_step_p_m5
                      * (
                              u(i - 1) - 2. * u(i) + u(i + 1)
                              + u_p(i - 1) - 2. * u_p(i) + u_p(i + 1)
                      )
                      * (
                              -u(i - 2) + 2. * u(i - 1) - 2. * u(i + 1) + u(i + 2)
                              - u_p(i - 2) + 2. * u_p(i - 1) - 2. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c001 / 16. * _x_step_p_m1
                      * (
                              u(i)
                              + u_p(i)
                      )
                      * (
                              u(i)
                              + u_p(i)
                      )
                      * (
                              u(i + 1) - u(i - 1)
                              + u_p(i + 1) - u_p(i - 1)
                      )
                    + c.c2 / 2. * _x_step_p_m2
                      * (u(i - 1) - 2. * u(i) + u(i + 1))
                    + c.c2 / 2. * _x_step_p_m2
                      * (u_p(i - 1) - 2. * u_p(i) + u_p(i + 1))
                    + c.c4 / 2. * _x_step_p_m4
                      * (
                              u(i - 2) - 4. * u(i - 1) + 6. * u(i) - 4. * u(i + 1) + u(i + 2)
                              + u_p(i - 2) - 4. * u_p(i - 1) + 6. * u_p(i) - 4. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c04 / 4. * _x_step_p_m4
                      * (
                              u(i)
                              + u_p(i)
                      )
                      * (
                              u(i - 2) - 4. * u(i - 1) + 6. * u(i) - 4. * u(i + 1) + u(i + 2)
                              + u_p(i - 2) - 4. * u_p(i - 1) + 6. * u_p(i) - 4. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c13 / 16. * _x_step_p_m4
                      * (
                              u(i + 1) - u(i - 1)
                              + u_p(i + 1) - u_p(i - 1)
                      )
                      * (
                              -u(i - 2) + 2. * u(i - 1) - 2. * u(i + 1) + u(i + 2)
                              - u_p(i - 2) + 2. * u_p(i - 1) - 2. * u_p(i + 1) + u_p(i + 2)
                      )
                    + c.c22 / 4. * _x_step_p_m4
                      * (
                              u(i - 1) - 2. * u(i) + u(i + 1)
                              + u_p(i - 1) - 2. * u_p(i) + u_p(i + 1)
                      )
                      * (
                              u(i - 1) - 2. * u(i) + u(i + 1)
                              + u_p(i - 1) - 2. * u_p(i) + u_p(i + 1)
                      );
    }

    system_rhs(n - 2) = (2. * u(n) - 5. * u(n - 1) + 4. * u(n - 2) - u(n - 3)) * _x_step_p_m2 - uxx_b(time);
    system_rhs(n - 1) = (3. * u(n) - 4. * u(n - 1) + u(n - 2)) / 2. * _x_step_p_m1 - ux_b(time);
    system_rhs(n) = u(n) - u_b(time);
}

std::ostream & KdV5Equation::dumpWithCoeffs(std::ostream &s, const KdV5Equation::coefficients &c) const {
    s << "    u_t";

    if (c.c3 > 0)
        s << " + " << c.c3 << " * u_xxx";
    else if (c.c3 < 0)
        s << " - " << -c.c3 << " * u_xxx";

    if (c.c5 > 0)
        s << " + " << c.c5 << " * u_xxxxx";
    else if (c.c5 < 0)
        s << " - " << -c.c5 << " * u_xxxxx";

    if (c.c01 > 0)
        s << " + " << c.c01 << " * u u_x";
    else if (c.c01 < 0)
        s << " - " << -c.c01 << " * u u_x";

    if (c.c03 > 0)
        s << " + " << c.c03 << " * u u_xxx";
    else if (c.c03 < 0)
        s << " - " << -c.c03 << " * u u_xxx";

    if (c.c12 > 0)
        s << " + " << c.c12 << " * u_x u_xx";
    else if (c.c12 < 0)
        s << " - " << -c.c12 << " * u_x u_xx";

    if (c.c05 > 0)
        s << " + " << c.c05 << " * u u_xxxxx";
    else if (c.c05 < 0)
        s << " - " << -c.c05 << " * u u_xxxxx";

    if (c.c14 > 0)
        s << " + " << c.c14 << " * u_x u_xxxx";
    else if (c.c14 < 0)
        s << " - " << -c.c14 << " * u_x u_xxxx";

    if (c.c23 > 0)
        s << " + " << c.c23 << " * u_xx u_xxx";
    else if (c.c23 < 0)
        s << " - " << -c.c23 << " * u_xx u_xxx";

    if (c.c001 > 0)
        s << " + " << c.c001 << " * u^2 u_x";
    else if (c.c001 < 0)
        s << " - " << -c.c001 << " * u^2 u_x";

    if (c.c2 > 0)
        s << " + " << c.c2 << " * u_xx";
    else if (c.c2 < 0)
        s << " - " << -c.c2 << " * u_xx";

    if (c.c4 > 0)
        s << " + " << c.c4 << " * u_xxxx";
    else if (c.c4 < 0)
        s << " - " << -c.c4 << " * u_xxxx";

    if (c.c04 > 0)
        s << " + " << c.c04 << " * u u_xxxx";
    else if (c.c04 < 0)
        s << " - " << -c.c04 << " * u u_xxxx";

    if (c.c13 > 0)
        s << " + " << c.c13 << " * u_x u_xxx";
    else if (c.c13 < 0)
        s << " - " << -c.c13 << " * u_x u_xxx";

    if (c.c22 > 0)
        s << " + " << c.c22 << " * u_xx^2";
    else if (c.c22 < 0)
        s << " - " << -c.c22 << " * u_xx^2";

    s << " = 0";

    return s;
}

std::ostream& operator<< (std::ostream& stream, const KdV5Equation& eq) {
    stream << "[equation]\n";
    eq.dumpWithCoeffs(stream, eq.coeffs);

    return stream;
}

std::map<double, double> KdV5Equation::computeKineticEnergy(unsigned int order, const expParameters &params) {
    std::map<double, double> energy_list;
    for (auto it = solution_list.begin(); it != solution_list.end(); it++) {
        VectorXd solution = it->second;
        const auto h = mesh.x_step;

        double energy;
        if (order == 0) {
            energy = solution.array().square().sum() * h;
        } else if (order == 1) {
            auto n = solution.size();
            energy = (
                             ((solution.tail(n - 2) - solution.head(n - 2)) / 2. / h).array().square() *
                             (params.beta - params.delta)
                             + ((solution.head(n - 2) - 2. * solution.segment(1, n - 2) + solution.tail(n - 2)) / h /
                                h).array().square() * params.gamma
                             + solution.array().cube() * params.epsilon / 2.
                     ).sum() * h;
        } else {
            throw std::logic_error("Order must be 0 or 1: " + order);
        }
        energy_list.insert(std::pair<double, double>(it->first, energy));
    }

    return energy_list;
}

std::map<double, double> KdV5Equation::computePotentialEnergy() {
    std::map<double, double> energy_list;
    for (auto it = solution_list.begin(); it != solution_list.end(); it++) {
        ArrayXd solution = (it->second).array();
        ArrayXd sign = solution / solution.abs();
        double energy = (solution.square() * sign).sum() * mesh.x_step;
        energy_list.insert(std::pair<double, double>(it->first, energy));
    }

    return energy_list;
}
