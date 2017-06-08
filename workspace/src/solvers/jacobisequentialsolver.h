#ifndef JACOBI_JACOBISEQUENTIALSOLVER_H
#define JACOBI_JACOBISEQUENTIALSOLVER_H

#include "jacobisolver.h"

class JacobiSequentialSolver : public JacobiSolver
{
public:
    JacobiSequentialSolver(const float** A, const float* b, int N) : JacobiSolver(A, b, N) {}

private:
    void deltax(const float *x, float *dest);
};

#endif //JACOBI_JACOBISEQUENTIALSOLVER_H
