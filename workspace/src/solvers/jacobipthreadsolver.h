//
// Created by caos on 06/06/17.
//

#ifndef JACOBI_JACOBIPTHREADSOLVER_H
#define JACOBI_JACOBIPTHREADSOLVER_H

#include <thread>
#include <queue>
#include "jacobisolver.h"

class JacobiPThreadSolver : public JacobiSolver
{
public:
    JacobiPThreadSolver(const double** A, const double* b, int N, int nWorkers);

private:
    int mnWorkers;

    void deltax(const double *x, double *dest);
    int getNWorkers() { return mnWorkers; };
    void threadBody(int minI, int maxI, const double* x, double* dest);

};


#endif //JACOBI_JACOBIPTHREADSOLVER_H
