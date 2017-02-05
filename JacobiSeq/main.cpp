#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace std;

#define MAX 50

#define VERBOSE
#define SHOW_TIME
#define READ_FROM_FILE

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

void vectDiff(const double* A, const double* B, double* C, int N)
{
    for(int i=0; i < N; i++)
        C[i] = A[i]-B[i];
}

void vectSum(const double* A, const double* B, double* C, int N)
{
    for(int i=0; i < N; i++)
        C[i] = A[i]+B[i];
}

void matvectProd(const double** A, const double* B, double* C, int N)
{
    for(int i=0; i < N; i++)
    {
        C[i] = 0.f;

        for (int j = 0; j < N; j ++)
            C[i] += A[i][j]*B[j];
    }
}

double sqnorm2(const double* A, int N)
{
    double sum = 0.f;
    for(int i = 0; i < N; i++)
        sum += A[i]*A[i];
    return sum;
}

/* Functions for Jacobi method */
void deltax(const double** D1, const double** A, const double* x, const double* b, double* res, int N)
{
    double* tmp = new double[N];

    matvectProd(A, x, tmp, N);
    vectDiff(b, tmp, tmp, N);
    matvectProd(D1, tmp, res, N);

    delete[] tmp;
}

bool converges(const double* deltax, int k, int N, int MAX_ITERATIONS, double EPS)
{
    return (k >= MAX_ITERATIONS || sqnorm2(deltax, N) <= EPS);
}

int main(int argc, char* argv[])
{
    /* Setup */
    srand(42);

    #ifdef READ_FROM_FILE
    int N, MAX_ITERATIONS;
    double EPS;
    ifstream in("input.in");

    in >> N;
    in >> MAX_ITERATIONS;
    in >> EPS;
    #else
    /* Read parameters */
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " N" << endl;
        return -1;
    }

    int N = atoi(argv[1]);
    int MAX_ITERATIONS = atoi(argv[2]);
    double EPS = atof(argv[3]);
    #endif

    /* Generate a diagonal dominant matrix of the correct size, and the vectors of known terms and solution */
    double** A = new double*[N];
    double* b = new double[N];
    double* x = new double[N];

    for (int i=0; i < N; i++)
        A[i] = new double[N];

    #ifdef READ_FROM_FILE
    for (int i=0; i < N; i++)
        for (int j = 0; j < N; j++)
            in >> A[i][j];

    for (int i=0; i < N; i++)
        in >> b[i];

    for (int i=0; i < N; i++)
        in >> x[i];

    in.close();
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

    /* Compute D^-1 and first deltax*/
    double** D1 = new double*[N];
    double* dx = new double[N];

    for (int i=0; i < N; i++)
    {
        D1[i] = new double[N];
        D1[i][i] = 1/A[i][i];
    }

    deltax((const double**)D1, (const double**)A, x, b, dx, N);

    /* Jacobi Algorithm */
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    int k = 0;
    while (!converges((const double*)dx, k, N, MAX_ITERATIONS, EPS))
    {
        deltax((const double**)D1, (const double**)A, x, b, dx, N);
        vectSum((const double*)x, (const double*)dx, x, N);
        k++;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedTime = end-start;

    #ifdef VERBOSE
    cout << endl << "## (x_solution^T) : ";
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";

    cout << endl << "2-norm error: " << sqnorm2(dx, N) << endl;
    #endif

    #ifdef SHOW_TIME
    cout << "Solved in " << k << " iterations." << endl;
    cout << "Took " << elapsedTime.count() << " s." << endl;
    #endif

    /* Release resources */
    for (int i=0; i < N; i++)
        delete[] A[i];

    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}