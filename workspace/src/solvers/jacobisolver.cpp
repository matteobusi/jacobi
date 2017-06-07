//
// Created by caos on 07/02/17.
//

#include <chrono>
#include <cstring>

#include "jacobisolver.h"

JacobiReport JacobiSolver::solve(int maxIterations, double eps, double *x)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;  
    JacobiReport jr;

    double* dx = new double[mN];

    int k = 0;
    bool conv = false;

    jr.convTime = 0.f;
    jr.compTime = 0.f;
    jr.updateTime = 0.f;

    while (!conv)
    {
        start = std::chrono::system_clock::now();
        deltax((const double*)x, dx);
        end = std::chrono::system_clock::now();

        jr.compTime += (end-start).count();

        start = std::chrono::system_clock::now();
        update(x, (const double*)dx);
        end = std::chrono::system_clock::now();

        jr.updateTime += (end-start).count();

        start = std::chrono::system_clock::now();
        conv = convergenceCheck(k, maxIterations, eps, (const double*)dx);
        end = std::chrono::system_clock::now();

        jr.convTime += (end-start).count();
        k++;
    }

    jr.error = norm(dx);
    jr.nIterations = k;
    jr.nWorkers = getNWorkers();

    delete[] dx;

    return jr;
}

void JacobiSolver::update(double *x, const double *dx)
{
    memcpy(x, dx, sizeof(double)*mN);
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