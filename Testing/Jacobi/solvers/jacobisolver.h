#ifndef JACOBI_JACOBI_H
#define JACOBI_JACOBI_H

#include "jacobireport.h"

class JacobiSolver
{
protected:
    const double** mA;
    const double* mb;
    int mN;

    virtual void deltax(const double* x, double *dest) = 0;
    virtual void update(double* x, const double* dx);
    virtual bool convergenceCheck(int k, int maxIterations, double eps, const double* dx);
    virtual double norm(const double *dx);

    virtual int getNWorkers() { return 1; }

public:
    JacobiSolver(const double** A, const double* b, int N) : mA(A), mb(b), mN(N) {};

    JacobiReport solve(int maxIterations, double eps, double* x);
};


#endif //JACOBI_JACOBI_H
