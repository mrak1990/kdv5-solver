#include "kdv_5.h"

KdV5Equation::meshProps::meshProps(double x_min, double x_max, double t_max, unsigned int mesh_size) :
        x_min(x_min), x_max(x_max), t_min(0.), t_max(t_max),
        x_step((x_max - x_min) / mesh_size), t_step(x_step),
        x_size(mesh_size), t_size(static_cast<unsigned int>((t_max - t_min) / t_step)) {}

KdV5Equation::meshProps::meshProps(double x_min, double x_max, double t_max, int mesh_size, double step_ratio) :
        x_min(x_min), x_max(x_max), t_min(0.), t_max(t_max),
        x_step((x_max - x_min) / mesh_size), t_step(step_ratio * x_step),
        x_size(mesh_size), t_size(static_cast<unsigned int>((t_max - t_min) / t_step)) {}

std::ostream& operator<< (std::ostream& stream, const KdV5Equation::meshProps& props) {
    stream << "    x min:  " << props.x_min << "\n"
           << "    x max:  " << props.x_max << "\n"
           << "    t min:  " << props.t_min << "\n"
           << "    t max:  " << props.t_max << "\n"
           << "    x step: " << props.x_step << "\n"
           << "    t step: " << props.t_step << "\n"
           << "    x size: " << props.x_size << "\n"
           << "    t size: " << props.t_size << std::endl;

    return stream;
}