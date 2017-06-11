#include "jacobisolver.h"

JacobiReport JacobiSolver::solve(int maxIterations, float eps, float *x)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;  
    std::chrono::time_point<std::chrono::system_clock> startL, endL;  

    int k = 0;
    bool conv = false;

    jr.convTime = 0.f;
    jr.updateTime = 0.f;

    // This is updated directly by deltax method!
    jr.sbTime = 0.f;

    startL = std::chrono::system_clock::now();
    while (!conv)
    {
        deltax((const float*)x, dx);

        start = std::chrono::system_clock::now();
        update(x, (const float*) dx);
        end = std::chrono::system_clock::now();

        jr.updateTime += ((std::chrono::duration<float>)(end-start)).count();
        k++;

        start = std::chrono::system_clock::now();
        conv = convergenceCheck(k, maxIterations, eps, (const float*)dx);
        end = std::chrono::system_clock::now();

        jr.convTime += ((std::chrono::duration<float>)(end-start)).count();
    }
    endL = std::chrono::system_clock::now();


    jr.error = norm(dx);
    jr.nIterations = k;
    jr.nWorkers = getNWorkers();
    jr.latency = ((std::chrono::duration<float>)(endL-startL)).count();

    return jr;
}

void JacobiSolver::update(float *x, const float *dx)
{
    memcpy(x, dx, sizeof(float)*mN);
}

bool JacobiSolver::convergenceCheck(int k, int maxIterations, float eps, const float *dx)
{
    return (k >= maxIterations || norm(dx) <= eps);
}

float JacobiSolver::norm(const float *dx)
{
    float norm = 0.f;
    for (int i=0; i < mN; i++)
        norm += dx[i]*dx[i];// fabs(dx[i]);

    return norm;
}