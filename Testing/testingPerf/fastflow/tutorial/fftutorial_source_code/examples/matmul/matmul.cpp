/* -*- Mode: C++; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ***************************************************************************
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2 as 
 *  published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  As a special exception, you may use this file as part of a free software
 *  library without restriction.  Specifically, if other files instantiate
 *  templates or use macros or inline functions from this file, or you compile
 *  this file and link it with other files to produce an executable, this
 *  file does not by itself cause the resulting executable to be covered by
 *  the GNU General Public License.  This exception does not however
 *  invalidate any other reasons why the executable file might be covered by
 *  the GNU General Public License.
 *
 ****************************************************************************
 */

/* 
 * Author: Massimo Torquati <torquati@di.unipi.it> 
 * Date:   July 2014
 */

/*
 * NxM matrix multiplication using the standard ijk algorithm O(n^3).
 *
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h>
#include <sys/time.h>

#include <ff/parallel_for.hpp>
using namespace ff;

const double THRESHOLD = 0.001;
bool  check    = false;  // set to true for checking the result
bool PFSPINWAIT=true;    // enabling spinWait
int  PFWORKERS =1;       // parallel_for parallelism degree
int  PFGRAIN   =0;       // default static scheduling of iterations


void random_init (long M, long N, long P, double *A, double *B) {
    for (long i = 0; i < M; i++)
        for (long j = 0; j < P; j++) {  
            A[i*P+j] = 5.0 - ((double)(rand()%100) / 10.0); 
        }
    for (long i = 0; i < P; i++) 
        for (long j = 0; j < N; j++) {  
            B[i*N+j] = 5.0 - ((double)(rand()%100) / 10.0);
        }     
}
    
long parse_arg(const char *string) {
    char *endptr;
    long r = strtol(string, &endptr, 10);
    assert(endptr[0]==0); // used up the whole string
    assert(r>0);          // need a positive integer
    return r;
}

void printarray(const double *A, long m, long n, long N) {
    for (long i=0; i<m; i++) {
        for (long j=0; j<n; j++)
            printf("%f\t", A[i*N+j]);
        printf("\n");
    }
}

long CheckResults(long m, long n, const double *C, const double *C1) {
	for (long i=0; i<m; i++)
		for (long j=0; j<n; j++) {
			long idx = i*n+j;
			if (fabs(C[idx] - C1[idx]) > THRESHOLD) {
				printf("ERROR %ld,%ld %f != %f\n", i, j, C[idx], C1[idx]);
				return 1;
			}
		}
	printf("OK.\n");
	return 0;
}

void get_time( void (*F)(long m, long n, long p, const double* A, long AN, const double* B, long BN, double* C, long CN),
               long m, long n, long p, const double* A, long AN, const double* B, long BN, double *C,
               const char *descr) {
	printf("Executing %-40s", descr);
	fflush(stdout);
	struct timeval before,after;
	gettimeofday(&before, NULL);
    
	F(m, n, p, A, AN, B, BN, C, BN);
	gettimeofday(&after,  NULL);
    
	double tdiff = after.tv_sec - before.tv_sec + (1e-6)*(after.tv_usec - before.tv_usec);
	printf("  Done in %11.6f secs.\n", tdiff);
}

void get_time_and_check(void (*F)(long m,long n,long p, const double* A, long AN, const double* B, long BN, double* C, long CN),
                        long m,long n,long p, const double* A, long AN, const double* B, long BN, const double *C_expected, double *C,
                        const char *descr) {
    get_time(F, m, n, p, A, AN, B, BN, C, descr);
    CheckResults(m, n, C_expected, C);
}

// triple nested loop (ijk) implementation 
void seqMatMult(long m, long n, long p, 
                const double* A, const long AN, 
                const double* B, const long BN, 
                double* C, const long CN)  {   

    for (long i = 0; i < m; i++)  
        for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;  
            for (long k = 0; k < p; k++)  
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
        }  
} 

void PFMatMultK(long m, long n, long p, 
                const double* A, const long AN, 
                const double* B, const long BN, 
                double* C, const long CN)  {   

    ParallelForReduce<double> pf(PFWORKERS, PFSPINWAIT);    
    auto Fsum = [](double& v, const double elem) { v += elem; };        
    for (long i = 0; i < m; i++)  
        for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;  
            pf.parallel_reduce(C[i*CN+j], 0.0,
                               0,p,[A,B,i,j,CN,AN,BN,&C](const long k, double &c) { //for (long k = 0; k < p; k++) {
                c += A[i*AN+k]*B[k*BN+j];  
                               },Fsum);
        }  
} 


void PFMatMultJ(long m, long n, long p, 
                const double* A, const long AN, 
                const double* B, const long BN, 
                double* C, const long CN)  {   

    ParallelFor pf(PFWORKERS, PFSPINWAIT);    
    for (long i = 0; i < m; i++)  
        pf.parallel_for(0,n, [A,B,i,CN,AN,BN,p,&C](const long j) { //for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;  
            for (long k = 0; k < p; k++) {
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
            }
            });
} 

void PFMatMultI(long m, long n, long p, 
                const double* A, const long AN, 
                const double* B, const long BN, 
                double* C, const long CN)  {   

    ParallelFor pf(PFWORKERS, PFSPINWAIT);
    pf.parallel_for(0,m,[A,B,CN,AN,BN,p,n,&C](const long i) { // for (long i = 0; i < m; i++) {
        for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;  
            for (long k = 0; k < p; k++) {
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
            }
        }
        });
} 


int main(int argc, char* argv[]) {     
    if (argc < 4) {
        printf("\n\tuse: %s M N P pfworkers:chunksize [check=0]\n", argv[0]);
        printf("\t       A is M by P\n");
        printf("\t       B is P by N\n");
        printf("\t check!=0 executes also the sequential ijk loops for checking the result\n\n");
        return -1;
    }

    long M = parse_arg(argv[1]);
    long N = parse_arg(argv[2]);
    long P = parse_arg(argv[3]);
    if (argc >= 5) {
        std::string pfarg(argv[4]);
        int n = pfarg.find_first_of(":");
        if (n>0) {
            PFWORKERS = atoi(pfarg.substr(0,n).c_str());
            PFGRAIN   = atoi(pfarg.substr(n+1).c_str());
        } else PFWORKERS = atoi(argv[5]);
    }
    if (argc >= 6) check = (atoi(argv[5])?true:false);

    const double *A = (double*)malloc(M*P*sizeof(double));
    const double *B = (double*)malloc(P*N*sizeof(double));
    assert(A); assert(B);

    random_init(M, N, P, const_cast<double*>(A), const_cast<double*>(B));
    double *C       = (double*)malloc(M*N*sizeof(double));
    
    if (check) {
        double *C2   = (double*)malloc(M*N*sizeof(double));
        get_time(seqMatMult, M,N,P, A, P, B, N, C, "Standard ijk algorithm");
        get_time_and_check(PFMatMultK, M, N, P, A, P, B, N, C, C2, "ParallelFor over K");
        get_time_and_check(PFMatMultJ, M, N, P, A, P, B, N, C, C2, "ParallelFor over J");
        get_time_and_check(PFMatMultI, M, N, P, A, P, B, N, C, C2, "ParallelFor over I");
        free(C2);
    } else {
        get_time(PFMatMultK, M, N, P, A, P, B, N, C, "ParallelFor over K");
        get_time(PFMatMultJ, M, N, P, A, P, B, N, C, "ParallelFor over J");
        get_time(PFMatMultI, M, N, P, A, P, B, N, C, "ParallelFor over I");
    }

    free((void*)A); free((void*)B); free(C);    
    return 0;  
} 
