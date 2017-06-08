#ifndef JACOBI_JACOBITHREADSOLVER_H
#define JACOBI_JACOBITHREADSOLVER_H

#include <thread>
#include <queue>
#include "jacobisolver.h"

class JacobiThreadSolver : public JacobiSolver
{
public:
    JacobiThreadSolver(const float **A, const float *b, int N, int nWorkers) : JacobiSolver(A, b, N), mnWorkers(nWorkers) {}

private:
    int mnWorkers;

    void deltax(const float *x, float *dest);
    int getNWorkers() { return mnWorkers; };
    void threadBody(int minI, int maxI, const float* x, float* dest);
};


#endif //JACOBI_JACOBITHREADSOLVER_H
