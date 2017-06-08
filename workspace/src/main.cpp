#include <iostream>
#include "solvers/jacobisolver.h"
#include "solvers/jacobisequentialsolver.h"
#include "solvers/jacobiffsolver.h"
#include "solvers/jacobithreadsolver.h"

using namespace std;

#define MAX 100

enum COMPUTATION_METHOD
{
    SEQUENTIAL,
    FASTFLOW,
    THREAD
};

/* Helper functions */
float randFloat(float min, float max)
{
    float v = ((float) rand()) / ((float) RAND_MAX);

    return v*(max-min) + min;
}

void fillDiagDominantMatrix(float** A, int r, int c)
{
    for (int i=0; i < r; i++)
    {
        float sum = 0;
        for (int j = 0; j < c; j++)
            if (j!=i)
            {
                float val = randFloat(-MAX, MAX);
                sum += abs(val);
                A[i][j] = val;
            }

        /* Change back A[i][i] to be > then sum(A[i][j]) */
        A[i][i] = sum + randFloat(1, MAX+1);
    }
}

void fillKnownVector(float* b, int r)
{
    for (int i=0; i < r; i++)
        b[i] = randFloat(-MAX, MAX);
}

void fillInitialGuessVector(float* x, int r)
{
    for (int i=0; i < r; i++)
        x[i] = randFloat(-MAX, MAX);
}

void usage(char* exec)
{
    cerr << "Usage: " << exec << " N ITER ERR METHOD [NWORKERS] [GRAIN]" << endl;
    cerr << "Where: " << endl;
    cerr << "N : is the size of the matrix A" << endl;
    cerr << "ITER : is the maximum number of iterations" << endl;
    cerr << "ERR : is the maximum norm of an acceptable error" << endl;
    cerr << "METHOD: is either" << endl;
    cerr << "\ts : indicating that the sequential implementation must be used" << endl;
    cerr << "\tf : indicating that the FastFlow implementation must be used" << endl;
    cerr << "\tt : indicating that the Thread implementation must be used" << endl;
    cerr << "NWORKERS : the number of workers that should be used (ignored if METHOD is 's')" << endl;
    cerr << "GRAIN : the grain of the computation (only if METHOD is 'f')" << endl;
    cerr << endl << "Produces a CSV line, in the form:" << endl;
    cerr << "\tN_WORKERS N_ITERATIONS COMP_TIME UPD_TIME CONV_TIME LATENCY ERROR" << endl;
}

int main(int argc, char* argv[])
{
    /* Setup */
    srand(42);

    /* Read parameters */
    if (argc < 5)
    {
        usage(argv[0]);
        exit(-1);
    }

    int N = atoi(argv[1]);
    int MAX_ITERATIONS = atoi(argv[2]);
    float EPS = atof(argv[3]);

    COMPUTATION_METHOD method = SEQUENTIAL;
    switch (argv[4][0]) /* First char is enough */
    {
        case 'f':
            method = FASTFLOW;
            break;
        case 't':
            method = THREAD;
            break;
    }
    
    int nworkers = 1;
    int grain = 1;

    if (method == THREAD)
        if (argc < 6)
        {
            usage(argv[0]);
            exit(-1);
        }
        else
            nworkers = atoi(argv[5]);
    else if (method == FASTFLOW)
        if (argc < 7)
        {
            usage(argv[0]);
            exit(-1);
        }
        else
        {
            nworkers = atoi(argv[5]);
            grain = atoi(argv[6]);
        }

    /* Generate a diagonal dominant matrix of the correct size, and the vectors of known terms and solution */
    float** A = new float*[N];
    float* b = new float[N];
    float* x = new float[N];

    if (A == NULL || b == NULL || x == NULL)
    {
        cerr << "Error while allocating resources." << endl;
        exit(-1);
    }

    for (int i=0; i < N; i++)
    {
        A[i] = new float[N];
        if (A[i] == NULL)
        {
            cerr << "Error while allocating resources." << endl;
            exit(-1);
        }
    }

    /* And fill them properly to assure convergence */
    fillDiagDominantMatrix(A, N, N);
    fillKnownVector(b, N);
    fillInitialGuessVector(x, N);

    /* Choose the proper solving method */
    JacobiSolver* js = NULL;

    if (method == SEQUENTIAL)
        js = new JacobiSequentialSolver((const float**)A, (const float*)b, N);
    else if (method == FASTFLOW)
        js = new JacobiFFSolver((const float**)A, (const float*)b, N, nworkers, grain);
    else if (method == THREAD)
        js = new JacobiThreadSolver((const float**)A, (const float*)b, N, nworkers);

    /* Solve the problem */
    JacobiReport report = js->solve(MAX_ITERATIONS, EPS, x);
    
    /* Produce results */
    cout << report << endl;

    /* Release resources */
    for (int i=0; i < N; i++)
        delete[] A[i];

    delete[] A;
    delete[] b;
    delete[] x;

    delete js;

    return 0;
}