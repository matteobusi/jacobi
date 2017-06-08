#include "jacobipthreadsolver.h"

void JacobiPThreadSolver::deltax(const float* x, float *dest)
{
    std::thread* workers = new std::thread[mnWorkers];

    int chunkSize = mN/mnWorkers;

    for (int w=0; w < mnWorkers; w++)
        workers[w] = std::thread(&JacobiPThreadSolver::threadBody, this, w*chunkSize, (w+1)*chunkSize, x, dest);

    for (int w=0; w < mnWorkers; w++)
        workers[w].join();

    delete [] workers;
}

void JacobiPThreadSolver::threadBody(int minI, int maxI, const float* x, float* dest)
{
    for (int i=minI; i < maxI; i++)
    {
        float s = 0.f;
        for (int j=0; j < i; j++)
            s+= mA[i][j]*x[j];

        for (int j=i+1; j < mN; j++)
            s+= mA[i][j]*x[j];
        dest[i] = (mb[i] - s)/mA[i][i];
    }
}

JacobiPThreadSolver::JacobiPThreadSolver(const float **A, const float *b, int N, int nWorkers) : JacobiSolver(A, b, N), mnWorkers(nWorkers)
{
}