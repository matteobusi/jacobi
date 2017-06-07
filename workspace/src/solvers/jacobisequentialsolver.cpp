#include "jacobisequentialsolver.h"

void JacobiSequentialSolver::deltax(const double* x, double *dest)
{
    // D^-1 ( b - Ax )
    for (int i=0; i < mN; i++)
    {
        double s = 0.f;
        for (int j=0; j < i; j++)
            s+= mA[i][j]*x[j];

        for (int j=i+1; j < mN; j++)
            s+= mA[i][j]*x[j];


        dest[i] = (mb[i] - s)/mA[i][i];
    }
}
