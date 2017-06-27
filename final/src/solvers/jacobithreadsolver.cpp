#include "jacobithreadsolver.h"

void JacobiThreadSolver::deltax(const float* x, float *dest)
{
    std::thread* workers = new std::thread[mnWorkers];
    int chunkSize = mN/mnWorkers;
    for (int w=0; w < mnWorkers; w++)
        workers[w] = std::thread(&JacobiThreadSolver::threadBody, this, w*chunkSize, (w+1)*chunkSize, x, dest);

    for (int w=0; w < mnWorkers; w++)
        workers[w].join();

    delete [] workers;
}

void JacobiThreadSolver::threadBody(int minI, int maxI, const float* x, float* dest)
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