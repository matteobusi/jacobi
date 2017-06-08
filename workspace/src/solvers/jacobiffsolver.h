//
// Created by caos on 07/02/17.
//

#ifndef JACOBI_JACOBIFFSOLVER_H
#define JACOBI_JACOBIFFSOLVER_H

#include <cmath>

#include <ff/parallel_for.hpp>

#include "jacobisolver.h"

class JacobiFFSolver : public JacobiSolver
{
public:
    JacobiFFSolver(const double** A, const double* b, int N, int nWorkers);
    virtual ~JacobiFFSolver() {delete pf;}

private:
    int mnWorkers;
    ff::ParallelFor* pf;
    void deltax(const double *x, double *dest);
    int getNWorkers() { return mnWorkers; };
};
#endif //JACOBI_JACOBIFFSOLVER_H
