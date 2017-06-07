#ifndef JACOBI_JACOBISEQUENTIALSOLVER_H
#define JACOBI_JACOBISEQUENTIALSOLVER_H

#include "jacobisolver.h"

class JacobiSequentialSolver : public JacobiSolver
{
public:
    JacobiSequentialSolver(const double** A, const double* b, int N) : JacobiSolver(A, b, N) {}

private:
    void deltax(const double *x, double *dest);
};

#endif //JACOBI_JACOBISEQUENTIALSOLVER_H
