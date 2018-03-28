#include "kdv_5.h"

KdV5Equation::coefficients::coefficients() {
    c3 = c5 = c01 = c03 = c12 = c05 = c14 = c23 = c001 = c2 = c4 = c04 = c13 = c22 = 0.;
}

std::string KdV5Equation::coefficients::dump() const {
    std::ostringstream str;
    str << "[equation]\n";
    str << "    u_t";

    if (c3 > 0)
        str << " + " << c3 << " * u_xxx";
    else if (c3 < 0)
        str << " - " << -c3 << " * u_xxx";

    if (c5 > 0)
        str << " + " << c5 << " * u_xxxxx";
    else if (c5 < 0)
        str << " - " << -c5 << " * u_xxxxx";

    if (c01 > 0)
        str << " + " << c01 << " * u u_x";
    else if (c01 < 0)
        str << " - " << -c01 << " * u u_x";

    if (c03 > 0)
        str << " + " << c03 << " * u u_xxx";
    else if (c03 < 0)
        str << " - " << -c03 << " * u u_xxx";

    if (c12 > 0)
        str << " + " << c12 << " * u_x u_xx";
    else if (c12 < 0)
        str << " - " << -c12 << " * u_x u_xx";

    if (c05 > 0)
        str << " + " << c05 << " * u u_xxxxx";
    else if (c05 < 0)
        str << " - " << -c05 << " * u u_xxxxx";

    if (c14 > 0)
        str << " + " << c14 << " * u_x u_xxxx";
    else if (c14 < 0)
        str << " - " << -c14 << " * u_x u_xxxx";

    if (c23 > 0)
        str << " + " << c23 << " * u_xx u_xxx";
    else if (c23 < 0)
        str << " - " << -c23 << " * u_xx u_xxx";

    if (c001 > 0)
        str << " + " << c001 << " * u^2 u_x";
    else if (c001 < 0)
        str << " - " << -c001 << " * u^2 u_x";


    if (c2 > 0)
        str << " + " << c2 << " * u_xx";
    else if (c2 < 0)
        str << " - " << -c2 << " * u_xx";

    if (c4 > 0)
        str << " + " << c4 << " * u_xxxx";
    else if (c4 < 0)
        str << " - " << -c4 << " * u_xxxx";

    if (c04 > 0)
        str << " + " << c04 << " * u u_xxxx";
    else if (c04 < 0)
        str << " - " << -c04 << " * u u_xxxx";

    if (c13 > 0)
        str << " + " << c13 << " * u_x u_xxx";
    else if (c13 < 0)
        str << " - " << -c13 << " * u_x u_xxx";

    if (c22 > 0)
        str << " + " << c22 << " * u_xx^2";
    else if (c22 < 0)
        str << " - " << -c22 << " * u_xx^2";

    return str.str();
}

std::ostream& operator<< (std::ostream& stream, const KdV5Equation::coefficients& params) {
    stream << "    c3=" << params.c3 << "\n"
           << "    c5=" << params.c5 << "\n"
           << "    c01=" << params.c01 << "\n"
           << "    c03=" << params.c03 << "\n"
           << "    c12=" << params.c12 << "\n"
           << "    c05=" << params.c05 << "\n"
           << "    c14=" << params.c14 << "\n"
           << "    c23=" << params.c23 << "\n"
           << "    c001=" << params.c001 << "\n"
           << "    c4=" << params.c4 << "\n"
           << "    c04=" << params.c04 << "\n"
           << "    c13=" << params.c13 << "\n"
           << "    c22=" << params.c22 << std::endl;

    return stream;
};