#include <iostream>

#include "solvers/jacobisolver.h"
#include "solvers/jacobisequentialsolver.h"
#include "solvers/jacobiffsolver.h"
#include "solvers/jacobipthreadsolver.h"

using namespace std;

#define MAX 50

// #define VERBOSE
#define PRINT_CSV
// #define READ_FROM_CIN

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

int main(int argc, char* argv[])
{
    /* Setup */
    srand(42);

    #ifdef READ_FROM_CIN
    int N, MAX_ITERATIONS;
    double EPS;
    int nworkers;
    COMPUTATION_METHOD method;

    cin >> N;
    cin >> MAX_ITERATIONS;
    cin >> EPS;
    cin >> method;
    cin >> nworkers;
    #else
    /* Read parameters */
    if (argc < 5)
    {
        cerr << "Usage: " << argv[0] << " N MAX_IT EPS METHOD [NWORKERS]" << endl;
        cerr << "Where " << "METHOD must be 0 for sequential, 1 for fast-flow and 2 for pthread." << endl;
        cerr << "Produces on stdout a line of a csv file, in the form: NWORKERS, NITERATIONS, TIME, ERROR" << endl;

        return -1;
    }

    int N = atoi(argv[1]);
    int MAX_ITERATIONS = atoi(argv[2]);
    double EPS = atof(argv[3]);
    COMPUTATION_METHOD method = (COMPUTATION_METHOD) atoi(argv[4]);
    int nworkers = 1;

    if (method != SEQUENTIAL)
        nworkers = atoi(argv[5]);
    #endif
    /* Generate a diagonal dominant matrix of the correct size, and the vectors of known terms and solution */
    double** A = new double*[N];
    double* b = new double[N];
    double* x = new double[N];

    for (int i=0; i < N; i++)
        A[i] = new double[N];

    #ifdef READ_FROM_CIN
    for (int i=0; i < N; i++)
        for (int j = 0; j < N; j++)
            cin >> A[i][j];

    for (int i=0; i < N; i++)
        cin >> b[i];

    for (int i=0; i < N; i++)
        cin >> x[i];
    #else
    fillDiagDominantMatrix(A, N, N);
    fillKnownVector(b, N);
    fillInitialGuessVector(x, N);
    #endif

    #ifdef VERBOSE
    cout << "## A ##" << endl;
    for (int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl << "## b^T : ";
    for (int i = 0; i < N; i++)
        cout << b[i] << " ";

    cout << endl << "## (x_0^T) : ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    #endif

    /* Jacobi Algorithm */
    JacobiSolver* js = NULL;

    if (method == SEQUENTIAL)
        js = new JacobiSequentialSolver((const double**)A, (const double*)b, N);
    else if (method == FASTFLOW)
        js = new JacobiFFSolver((const double**)A, (const double*)b, N, nworkers);
    else if (method == PTHREAD)
        js = new JacobiPThreadSolver((const double**)A, (const double*)b, N, nworkers);

    JacobiReport report = js->solve(MAX_ITERATIONS, EPS, x);

    #ifdef VERBOSE
    cout << endl << "## (x_solution^T) : ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";

    cout << endl << "Error: " << report.error << endl;
    cout << "Solved in " << report.nIterations << " iterations." << endl;
    cout << "Took " << report.time << " s using " << report.nWorkers << " workers." << endl;
    #endif

    #ifdef PRINT_CSV
    cout << report << endl;
    #endif

    /* Release resources */
    for (int i=0; i < N; i++)
        delete[] A[i];

    delete[] A;
    delete[] b;
    delete[] x;

    delete js;

    return 0;
}