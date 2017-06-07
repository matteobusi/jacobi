//
// Created by caos on 07/02/17.
//

#include <chrono>
#include <cstring>
#include "jacobisolver.h"

JacobiReport JacobiSolver::solve(int maxIterations, double eps, double *x)
{
    JacobiReport jr;

    double* dx = new double[mN];

    deltax((const double*)x, dx);

    int k = 0;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    while (!convergenceCheck(k, maxIterations, eps, (const double*)dx))
    {
        deltax((const double*)x, dx);
        update(x, (const double*)dx);
        k++;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime = end-start;

    jr.error = norm(dx);
    jr.nIterations = k;
    jr.nWorkers = getNWorkers();
    jr.time = elapsedTime.count();

    delete[] dx;

    return jr;
}

void JacobiSolver::update(double *x, const double *dx)
{
    memcpy(x, dx, mN * sizeof(double));
}

bool JacobiSolver::convergenceCheck(int k, int maxIterations, double eps, const double *dx)
{
    return (k >= maxIterations || norm(dx) <= eps);
}

double JacobiSolver::norm(const double *dx)
{
    double norm = 0.f;
    for (int i=0; i < mN; i++)
        norm += dx[i]*dx[i];// fabs(dx[i]);

    return norm;
}