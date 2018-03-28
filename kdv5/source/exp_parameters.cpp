#include "kdv_5.h"

KdV5Equation::expParameters::operator KdV5Equation::coefficients() const {
    KdV5Equation::coefficients coeffs;

    coeffs.c3 = 1. / 6. * mu - 1. / 2. * beta + 1. / 2. * delta;
    coeffs.c4 = 1. / 2. * tau;
    coeffs.c5 = -1. / 12. * mu * beta + 1. / 4. * mu * delta - 1. / 4. * beta * delta
                + 1. / 2. * gamma + 19. / 360. * mu * mu + 3. / 8. * delta * delta - 1. / 8. * beta * beta;
//    coeffs.c6 = 1./12.*mu*tau
    coeffs.c01 = 3. / 2. * epsilon;
    coeffs.c03 = 7. / 4. * epsilon * delta - 1. / 4. * epsilon * beta + 5. / 12. * mu * epsilon;
    coeffs.c12 = 5. / 8. * epsilon * beta + 23. / 24. * epsilon * mu + 31. / 8. * epsilon * delta;
    coeffs.c14 = 1. / 4. * epsilon * tau;
    coeffs.c22 = -9. / 8. * epsilon * tau;
    coeffs.c13 = -epsilon * tau;
    coeffs.c05 = 1. / 4. * epsilon * gamma;
    coeffs.c14 = -11. / 8. * epsilon * gamma;
    coeffs.c23 = -15. / 4. * epsilon * gamma;
    coeffs.c001 = -3. / 8. * epsilon * epsilon;
    coeffs.c04 = 1. / 4. * epsilon * tau;

    return coeffs;
}

std::ostream& operator<< (std::ostream& stream, const KdV5Equation::expParameters& params) {
    stream << "    epsilon=" << params.epsilon << "\n"
           << "    mu=" << params.mu << "\n"
           << "    gamma=" << params.gamma << "\n"
           << "    beta=" << params.beta << "\n"
           << "    delta=" << params.delta << "\n"
           << "    tau=" << params.tau << std::endl;

    return stream;
}
