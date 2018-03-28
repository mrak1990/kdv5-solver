#ifndef KDV5_SERVICE_H
#define KDV5_SERVICE_H

#include <fstream>
#include <map>
#include <tuple>
#include <functional>
#include <iostream>

#include <Eigen/Dense>

using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Infinity;

void writeSolutionToFile(const ArrayXd &mesh, std::map<double, VectorXd> solution_list);

void writeErrorToFile(const ArrayXd &x, const ArrayXd &y);

double getCurrentError(const ArrayXd &mesh, double time, const VectorXd &solution,
                       std::function<double(const double, const double)> exact_sol);

ArrayXd getErrorList(const ArrayXd &mesh, const std::map<double, VectorXd> &solution_list,
                     std::function<double(const double, const double)> exact_sol);

std::tuple<double, double> computeLinearFitting(const ArrayXd &x, const ArrayXd &y);

#endif //KDV5_SERVICE_H
