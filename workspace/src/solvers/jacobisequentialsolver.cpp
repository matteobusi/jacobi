#include "jacobisequentialsolver.h"

void JacobiSequentialSolver::deltax(const float* x, float *dest)
{
    // Just update the compTime! 
    std::chrono::time_point<std::chrono::system_clock> start, end;  

    start = std::chrono::system_clock::now();    
    for (int i=0; i < mN; i++)
    {
        float s = 0.f;
        for (int j=0; j < i; j++)
            s+= mA[i][j]*x[j];

        for (int j=i+1; j < mN; j++)
            s+= mA[i][j]*x[j];


        dest[i] = (mb[i] - s)/mA[i][i];
    }
    end = std::chrono::system_clock::now();

    jr.compTime += ((std::chrono::duration<float>)(end-start)).count();
}
