#ifndef JACOBI_JACOBIPTHREADSOLVER_H
#define JACOBI_JACOBIPTHREADSOLVER_H

#include <thread>
#include <queue>
#include "jacobisolver.h"

class JacobiPThreadSolver : public JacobiSolver
{
public:
    JacobiPThreadSolver(const float** A, const float* b, int N, int nWorkers);

private:
    int mnWorkers;

    void deltax(const float *x, float *dest);
    int getNWorkers() { return mnWorkers; };
    void threadBody(int minI, int maxI, const float* x, float* dest);

};


#endif //JACOBI_JACOBIPTHREADSOLVER_H
