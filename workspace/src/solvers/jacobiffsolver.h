#ifndef JACOBI_JACOBIFFSOLVER_H
#define JACOBI_JACOBIFFSOLVER_H

#include <cmath>

#include <ff/parallel_for.hpp>

#include "jacobisolver.h"

class JacobiFFSolver : public JacobiSolver
{
public:
    JacobiFFSolver(const float** A, const float* b, int N, int nWorkers, int grain);
    virtual ~JacobiFFSolver() {delete pf;}

private:
    int mnWorkers;
    int mGrain;
    ff::ParallelFor* pf;
    void deltax(const float *x, float *dest);
    int getNWorkers() { return mnWorkers; };
};
#endif //JACOBI_JACOBIFFSOLVER_H
