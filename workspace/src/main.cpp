#include <iostream>
#include "solvers/jacobisolver.h"
#include "solvers/jacobisequentialsolver.h"
#include "solvers/jacobiffsolver.h"
#include "solvers/jacobipthreadsolver.h"

using namespace std;

#define MAX 100

enum COMPUTATION_METHOD
{
    SEQUENTIAL,
    FASTFLOW,
    PTHREAD
};

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
        for (int j = 0; j < c; j++)
            if (j!=i)
            {
                double val = randDouble(-MAX, MAX);
                sum += abs(val);
                A[i][j] = val;
            }

        /* Change back A[i][i] to be > then sum(A[i][j]) */
        A[i][i] = sum + randDouble(1, MAX+1);
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

void usage(char* exec)
{
    cerr << "Usage: " << exec << " N ITER ERR METHOD [NWORKERS]" << endl;
    cerr << "Where: " << endl;
    cerr << "N : is the size of the matrix A" << endl;
    cerr << "ITER : is the maximum number of iterations" << endl;
    cerr << "ERR : is the maximum norm of an acceptable error" << endl;
    cerr << "METHOD: is either" << endl;
    cerr << "\ts : indicating that the sequential implementation must be used" << endl;
    cerr << "\tf : indicating that the FastFlow implementation must be used" << endl;
    cerr << "\tp : indicating that the PThread implementation must be used" << endl;
    cerr << "NWORKERS : the number of workers that should be used (ignored if METHOD is 's')" << endl;
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
    double EPS = atof(argv[3]);

    COMPUTATION_METHOD method = SEQUENTIAL;
    switch (argv[4][0]) /* First char is enough */
    {
        case 'f':
            method = FASTFLOW;
            break;
        case 'p':
            method = PTHREAD;
            break;
    }
    
    int nworkers = 1;

    if (method != SEQUENTIAL && argc > 5)
        nworkers = atoi(argv[5]);
    else if (method != SEQUENTIAL && argc < 6)
    {
        usage(argv[0]);
        exit(-1);
    }

    /* Generate a diagonal dominant matrix of the correct size, and the vectors of known terms and solution */
    double** A = new double*[N];
    double* b = new double[N];
    double* x = new double[N];

    if (A == NULL || b == NULL || x == NULL)
    {
        cerr << "Error while allocating resources." << endl;
        exit(-1);
    }

    for (int i=0; i < N; i++)
    {
        A[i] = new double[N];
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
        js = new JacobiSequentialSolver((const double**)A, (const double*)b, N);
    else if (method == FASTFLOW)
        js = new JacobiFFSolver((const double**)A, (const double*)b, N, nworkers);
    else if (method == PTHREAD)
        js = new JacobiPThreadSolver((const double**)A, (const double*)b, N, nworkers);

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