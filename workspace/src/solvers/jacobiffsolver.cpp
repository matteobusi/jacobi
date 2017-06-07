//
// Created by caos on 07/02/17.
//

#include "jacobiffsolver.h"

void JacobiFFSolver::deltax(const double* x, double *dest)
{
    // D^-1 ( b - Ax )
    ff::ParallelFor pf(mnWorkers, false);

    pf.parallel_for(0, mN, 1, 1, [&](const long i)
    {
        double s = 0.f;
        for (int j=0; j < i; j++)
            s+= mA[i][j]*x[j];

        for (int j=i+1; j < mN; j++)
            s+= mA[i][j]*x[j];
        dest[i] = (mb[i] - s)/mA[i][i];
    });
}