#include "kdv_5.h"

#include "lapacke.h"

KdV5Equation::DGBSV::DGBSV(int N) : N(N) {
    IPIV = new int[N];
}

KdV5Equation::DGBSV::~DGBSV() {
    delete[] IPIV;
}

void KdV5Equation::DGBSV::solve(double *AB, double *B) {
    const int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB);

    if (info != 0) {
        throw std::logic_error("dgbsv solver return non-zero INFO: " + info);
    }
}