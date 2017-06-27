#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include <chrono>
#include <cstring>

#include "jacobireport.h"

class JacobiSolver
{
protected:
    const float** mA;
    const float* mb;
    float* dx;
    int mN;
    JacobiReport jr;

    virtual void deltax(const float* x, float *dest) = 0;
    virtual void update(float* x, const float* dx);
    virtual bool convergenceCheck(int k, int maxIterations, float eps, const float* dx);
    virtual float norm(const float *dx);

    virtual int getNWorkers() { return 1; }

public:
    JacobiSolver(const float** A, const float* b, int N) : mA(A), mb(b), mN(N) { dx = new float[mN]; };
    virtual ~JacobiSolver() { delete[] dx; }
    JacobiReport solve(int maxIterations, float eps, float* x);
};


#endif //JACOBI_JACOBI_H
