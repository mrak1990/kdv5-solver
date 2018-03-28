#ifndef KDV5_KDV_5_H
#define KDV5_KDV_5_H

#include <iostream>
#include <map>

#include <Eigen/Dense>

using Eigen::ArrayXd;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class KdV5Equation {
public:
    typedef std::function<double(const double)> func;
    struct coefficients {
        double c3;      // u_xxx
        double c5;      // u_xxxxx
        double c01;     // u u_x
        double c03;     // u u_xxx
        double c12;     // u_x u_xx
        double c05;     // u u_xxxxx
        double c14;     // u_x u_xxxx
        double c23;     // u_xx u_xxx
        double c001;    // u^2 u_x

        double c2;   // u_xx
        double c4;      // u_xxxx
        double c04;     // u u_xxxx
        double c13;     // u_x u_xxx
        double c22;     // u_xx^2

        coefficients();

        std::string dump() const;
        friend std::ostream& operator<< (std::ostream& stream, const coefficients& params);
    };
    struct expParameters {
        double epsilon;
        double mu;
        double gamma;
        double beta;
        double delta;
        double tau;
        expParameters() = default;
        operator KdV5Equation::coefficients() const;
        friend std::ostream& operator<< (std::ostream& stream, const expParameters& params);
    };
    struct modelParameters {
        double a;
        double l;
        double h_0;
        double h_i;
        double sigma;
        double rho_i;
        double rho_w;
        double E;
        double nu;
        double g;
        modelParameters();
        operator KdV5Equation::expParameters() const;
        friend std::ostream& operator<< (std::ostream& stream, const modelParameters& params);
    };
    struct meshProps {
        meshProps(double x_min, double x_max, double t_max, unsigned int mesh_size);
        meshProps(double x_min, double x_max, double t_max, int mesh_size, double step_ratio);
        double x_min;
        double x_max;
        double t_min;
        double t_max;
        double x_step;
        double t_step;
        int x_size;
        int t_size;
        friend std::ostream& operator<< (std::ostream& stream, const meshProps& props);
    };
    struct iterationsProps {
        unsigned int max;
        double max_error;
        unsigned int exact;
    };
    struct outProps {
        std::ostream* debug_stream = nullptr;
        int timestep_per_plot;
        int timestep_number_to_plot;
        std::string solution_file;
        bool instant_dump = true;
    };

    KdV5Equation(coefficients coefficients,
                 meshProps mesh2, iterationsProps iterations, unsigned int plot_count,
                 func initial, func u_a, func ux_a,
                 func uxx_a, func u_b, func ux_b, func uxx_b);
    void run();
    void setDebugStream(std::ostream &stream);
    void setInstantDump(bool instant);
    void setSolutionFile(std::string file);
    const ArrayXd &getXGrid() const;
    friend std::ostream& operator<< (std::ostream& stream, const KdV5Equation& props);
    const std::map<double, VectorXd> &getSolutions() const;
    std::map<double, double> computeKineticEnergy(unsigned int order, const expParameters &params);
    std::map<double, double> computePotentialEnergy();
private:
    class DGBSV {
    public:
        explicit DGBSV(int N);
        ~DGBSV();
        void solve(double *AB, double *B);
    private:
        int *IPIV;
        const int N;
        const int KU = 3;
        const int KL = 3;
        const int NRHS = 1;
        const int LDAB = 2 * KL + KU + 1;
        const int LDB = N;
    };

private:
    void assemble_system();
    void assemble_rhs();
    void dumpSolution();
    std::ostream &dumpWithCoeffs(std::ostream &s, const KdV5Equation::coefficients &c) const;

    DGBSV solver;
    const coefficients coeffs;
    meshProps mesh;
    const iterationsProps iterations;
    outProps out;
    const unsigned int plot_count;
    const func initial, u_a, ux_a, uxx_a, u_b, ux_b, uxx_b;
    ArrayXd xGrid;
    VectorXd solution_u, old_solution_u;
    MatrixXd system_matrix;
    VectorXd system_rhs;
    double time;
    std::map<double, VectorXd> solution_list;
    double _x_step_p_m1, _x_step_p_m2, _x_step_p_m3, _x_step_p_m4, _x_step_p_m5, _t_step_p_m1;
};

#endif //KDV5_KDV_5_H
