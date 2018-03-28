#include "kdv_5.h"

KdV5Equation::modelParameters::modelParameters() {
    a = l = h_0 = h_i = 0.;
    sigma = 1E4;
    rho_i = 900;
    rho_w = 1025;
    E = 3E9;
    nu = 0.3;
    g = 9.8;
}

KdV5Equation::modelParameters::operator KdV5Equation::expParameters() const {
    if (a == 0. || l == 0. || h_0 == 0.) {
        throw std::logic_error("Model parameters must be initialized");
    }

    KdV5Equation::expParameters params = {};
    params.epsilon = a / h_0;
    params.mu = h_0 * h_0 / l / l;
    params.gamma = E * std::pow(h_i, 3) / 12. / (1 - nu * nu) / rho_w / g / std::pow(l, 4);
    params.beta = h_i * sigma / rho_w / g / l / l;
    params.delta = rho_i / rho_w * h_i * h_0 / l / l;
    params.tau = 3E-3;
//    params.tau = 0;

    return params;
}

std::ostream& operator<< (std::ostream& stream, const KdV5Equation::modelParameters& params) {
    stream << "    a=" << params.a << "\n"
           << "    l=" << params.l << "\n"
           << "    h_0=" << params.h_0 << "\n"
           << "    h_i=" << params.h_i << "\n"
           << "    sigma=" << params.sigma << "\n"
           << "    rho_i=" << params.rho_i << "\n"
           << "    rho_w=" << params.rho_w << "\n"
           << "    E=" << params.E << "\n"
           << "    nu=" << params.nu << "\n"
           << "    g=" << params.g << std::endl;

    return stream;
};
