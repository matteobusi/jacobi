#include <iostream>
#include <chrono>
#include <fstream>

#include <ff/parallel_for.hpp>

using namespace std;
using namespace ff;

#define MAX 50

// #define VERBOSE
#define SHOW_TIME
// #define READ_FROM_FILE

/* Helper functions */
double randDouble(double min, double max)
{
    double v = ((double) rand()) / ((double) RAND_MAX);

    return v*(max-min) + min;
}

void fillDiagDominantMatrix(double** A, int r, int c)
{
    for (int i=0; i < r; i++)
    {
        double sum = 0;
        for (int j = 0; j < c; j ++)
            if (j!=i)
            {
                double val = randDouble(0, MAX);
                sum += val;
                A[i][j] = val;
            }

        /* Change back A[i][i] to be > then sum(A[i][j]) */
        A[i][i] = randDouble(sum+1, MAX+sum+1);
    }
}

void fillKnownVector(double* b, int r)
{
    for (int i=0; i < r; i++)
        b[i] = randDouble(-MAX, MAX);
}

void fillInitialGuessVector(double* x, int r)
{
    for (int i=0; i < r; i++)
        x[i] = randDouble(-MAX, MAX);
}


void matvectProdP(const double** A, const double* B, double* C, int N, int nworkers)
{
    ParallelFor pf(nworkers, false);

    pf.parallel_for(0,N,[&](const long i) {
        C[i] = 0.f;

        for (int j = 0; j < N; j ++)
            C[i] += A[i][j]*B[j];
    });
}

void matvectProdNP(const double** A, const double* B, double* C, int N)
{
    for(int i=0; i < N; i++)
    {
        C[i] = 0.f;

        for (int j = 0; j < N; j ++)
            C[i] += A[i][j]*B[j];
    }
}

int main(int argc, char* argv[])
{
    /* Setup */
    srand(42);

    /* Read parameters */
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " N MAX_ITERATIONS nworkers" << endl;
        return -1;
    }

    int N = atoi(argv[1]);
    int MAX_ITERATIONS = atoi(argv[2]);
    int nworkers = atoi(argv[3]);

    /* Generate a diagonal dominant matrix of the correct size, and the vectors of known terms and solution */
    double** A = new double*[N];
    double* b = new double[N];
    double* x = new double[N];

    for (int i=0; i < N; i++)
        A[i] = new double[N];

    fillDiagDominantMatrix(A, N, N);
    fillKnownVector(b, N);
    fillInitialGuessVector(x, N);

    double sumP = 0.f;
    double sumNP = 0.f;

    for (int i=0; i < MAX_ITERATIONS; i++) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        matvectProdP((const double**)A, x, b, N, nworkers);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsedTime = end - start;
        sumP += elapsedTime.count();

        start = std::chrono::system_clock::now();

        matvectProdNP((const double**)A, x, b, N);

        end = std::chrono::system_clock::now();
        elapsedTime = end - start;
        sumNP += elapsedTime.count();
    }

    cout << "ParFor mean time: " << sumP/MAX_ITERATIONS << endl;
    cout << "Seq mean time: " << sumNP/MAX_ITERATIONS << endl;

    /* Release resources */
    for (int i=0; i < N; i++)
        delete[] A[i];

    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}