//
// Created by caos on 07/02/17.
//

#include "jacobipthreadsolver.h"

void JacobiPThreadSolver::deltax(const double* x, double *dest)
{
    std::thread* workers = new std::thread[mnWorkers];

    int chunkSize = mN/mnWorkers;

    for (int w=0; w < mnWorkers; w++)
        workers[w] = std::thread(&JacobiPThreadSolver::threadBody, this, w*chunkSize, (w+1)*chunkSize, x, dest);

    for (int w=0; w < mnWorkers; w++)
        workers[w].join();

    delete [] workers;
}

void JacobiPThreadSolver::threadBody(int minI, int maxI, const double* x, double* dest)
{
    for (int i=minI; i < maxI; i++)
    {
        double s = 0.f;
        for (int j=0; j < i; j++)
            s+= mA[i][j]*x[j];

        for (int j=i+1; j < mN; j++)
            s+= mA[i][j]*x[j];
        dest[i] = (mb[i] - s)/mA[i][i];
    }
}

JacobiPThreadSolver::JacobiPThreadSolver(const double **A, const double *b, int N, int nWorkers) : JacobiSolver(A, b, N), mnWorkers(nWorkers)
{
}