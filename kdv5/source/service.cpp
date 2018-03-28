#include "service.h"

void writeSolutionToFile(const std::string& file, const ArrayXd &mesh, std::map<double, VectorXd> solution_list) {
    std::ofstream out;
    out.open(file);

    for (auto it = solution_list.begin(); it != solution_list.end(); it++) {
        MatrixXd m(mesh.size(), 3);
        m.col(0).setConstant(it->first);
        m.col(1) = mesh;
        m.col(2) = it->second;

        out << m << std::endl << std::endl << std::endl;
    }

    out.close();
}

double getCurrentError(const ArrayXd &mesh, double time, const VectorXd &solution,
                       std::function<double(const double, const double)> exact_sol) {
    VectorXd exact(mesh.size());
    for (unsigned int i = 0; i <= mesh.size() - 1; i++)
        exact(i) = exact_sol(mesh(i), time);

    return (exact - solution).lpNorm<Infinity>();
}

ArrayXd getErrorList(const ArrayXd &mesh, const std::map<double, VectorXd> &solution_list,
                     std::function<double(const double, const double)> exact_sol) {
    ArrayXd errors(solution_list.size());
    int j = 0;
    for (auto it = solution_list.begin(); it != solution_list.end(); it++, j++) {
        errors(j) = getCurrentError(mesh, it->first, it->second, exact_sol);
    }

    return errors;
}

std::tuple<double, double> computeLinearFitting(const ArrayXd &x, const ArrayXd &y) {
    const ArrayXd xy = x * y;
    const ArrayXd x2 = x.square();
    const double x2_mean = x2.mean();
    const double xy_mean = xy.mean();
    const double x_mean = x.mean();
    const double y_mean = y.mean();

    const double a = (xy_mean - x_mean * y_mean) / (x2_mean - x_mean * x_mean);
    const double b = y_mean - a * x_mean;

    return std::make_tuple(a, b);
}