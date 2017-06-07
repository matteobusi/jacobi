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
 * Date:   September 2014
 */

/*
 * NxM matrix multiplication using the block-based algorithm
 * (the implementation is not optimised, it doesn't store the blocks
 * in a row-major format).
 *
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h>
#include <sys/time.h>
#include <ff/taskf.hpp>
using namespace ff;

const double THRESHOLD = 0.001;
bool  check = false;  // set to true for checking the result

void random_init (long M, long N, long P, double *A, double *B) {
    for (long i = 0; i < M; i++) {  
        for (long j = 0; j < P; j++) {  
            A[i*P+j] = 5.0 - ((double)(rand()%100) / 10.0); 
        }
    }
	
    for (long i = 0; i < P; i++) {  
        for (long j = 0; j < N; j++) {  
            B[i*N+j] = 5.0 - ((double)(rand()%100) / 10.0);
        }     
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


// triple nested loop (ijk) implementation 
void seqMatMult(ff_taskf *, long notused, 
                long m, long n, long p, 
                const double* A, const long AN, 
                const double* B, const long BN, 
                double* C, const long CN)  {   
    (void)notused;
    for (long i = 0; i < m; i++)  
        for (long j = 0; j < n; j++) {
            C[i*CN+j] = 0.0;  
            for (long k = 0; k < p; k++)  
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
        }  
} 

void blockMM(long b, const double *A, const long AN, 
             const double *B, const long BN, 
             double *C, const long CN) {

    for (long i = 0; i < b; i++)  
        for (long j = 0; j < b; j++) {
            for (long k = 0; k < b; k++)  
                C[i*CN+j] += A[i*AN+k]*B[k*BN+j];  
        }  
}

void zerosBlock(long b, double *C, const long CN) {
    for (long i = 0; i < b; i++)  
        for (long j = 0; j < b; j++) 
            C[i*CN+j] = 0.0;
}


/* 
 * Block-based algorithm:
 *  given a block size b and m = 2*b,  n = 3*b,  p = 4*b
 *
 *                            
 *  -----------------------     -----------------      -----------------------------
 * | A11 | A12 | A13 | A14 |   | B11 | B21 | B23 |    | A11*B11 | A11*B12 | A11*B13 |
 * |-----|-----|-----|-----| X |-----|-----|-----| =  | --------|---------|---------| +
 * | A21 | A22 | A23 | A24 |   | B21 | B22 | B23 |    | A21*B11 | A21*B12 | A21*B13 |
 *  -----------------------    |-----|-----|-----|     -----------------------------
 *                             | B31 | B32 | B33 |       -----------------------------
 *                             |-----|-----|-----|      | A12*B21 | A12*B22 | A12*B23 |
 *                             | B41 | B42 | B43 |      | --------|---------|---------| +
 *                              -----------------       | A22*B21 | A22*B22 | A22*B23 |
 *  Algo:                                                -----------------------------
 *                                                         -----------------------------
 *  - fill out the C matrix with zeros                    | A13*B31 | A13*B32 | A13*B33 |
 *  - for k 1->pblocks                                    | --------|---------|---------| +
 *      for i 1->mblocks                                  | A23*B31 | A23*B32 | A23*B33 |
 *        for j 1->nblocks                                 -----------------------------
 *           Ci,j += Ai,k * Bk,j                             -----------------------------
 *                                                          | A14*B41 | A14*B42 | A14*B43 |
 *                                                          | --------|---------|---------|
 *                                                          | A24*B41 | A24*B42 | A24*B43 |
 *                                                           -----------------------------
 */
void blockMMult (ff_taskf *taskf, long b, 
                 long m, long n, long p,
                 const double *A, long AN,  // input:  m by p with row stride AN
                 const double *B, long BN,  // input:  p by n with row stride BN
                 double *C, long CN) {      // output: m by n with row stride CN.

	long mblocks = m/b;
	long nblocks = n/b;    
	long pblocks = p/b;

    taskf->run(); // starts the task scheduler 
    for(long i=0; i<mblocks; ++i) 
        for(long j=0; j<nblocks; ++j) {
            auto C_ij = [CN,i,j,b, &C]() {
                double *Cij       = &C[b*(i*CN+j)];
                zerosBlock(b,Cij,CN);
            };
            taskf->AddTask(C_ij);
        }
    taskf->wait();
    taskf->run();
    for(long i=0; i<mblocks; ++i) {
        for(long j=0; j<nblocks; ++j) {
            auto C_ij = [AN,BN,CN,i,j,b,pblocks, &A,&B,&C]() {
                for(long k=0; k<pblocks; ++k) {
                    const double *Aik = &A[b*(i*AN + k)];
                    const double *Bkj = &B[b*(k*BN + j)];
                    double *Cij       = &C[b*(i*CN+j)];
                    blockMM(b, Aik, AN, Bkj, BN, Cij, CN);                        
                }
            };            
            taskf->AddTask(C_ij); // enqueue a task
        }
    }
    taskf->wait(); // barrier        
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

void get_time( void (*F)(ff_taskf *, long b, long m, long n, long p, const double* A, long AN, const double* B, long BN, double* C, long CN), ff_taskf *taskf,
               long b, long m, long n, long p, const double* A, long AN, const double* B, long BN, double *C,
               const char *descr) {
	printf("Executing %-40s", descr);
	fflush(stdout);
	struct timeval before,after;
	gettimeofday(&before, NULL);
    
	F(taskf, b, m, n, p, A, AN, B, BN, C, BN);
	gettimeofday(&after,  NULL);
    
	double tdiff = after.tv_sec - before.tv_sec + (1e-6)*(after.tv_usec - before.tv_usec);
	printf("  Done in %11.6f secs.\n", tdiff);
}

void get_time_and_check(void (*F)(ff_taskf*, long b, long m,long n,long p, const double* A, long AN, const double* B, long BN, double* C, long CN), ff_taskf *taskf,
                        long b, long m,long n,long p, const double* A, long AN, const double* B, long BN, const double *C_expected, double *C,
                        const char *descr) {
    get_time(F, taskf, b, m, n, p, A, AN, B, BN, C, descr);
    CheckResults(m, n, C_expected, C);
}


int main(int argc, char* argv[]) {     
    if (argc < 6) {
        printf("\n\tuse: %s <M> <N> <P> <blocksize> <nworkers> [check=0]\n", argv[0]);
        printf("\t       <-> required argument, [-] optional argument\n");
        printf("\t       A is M by P\n");
        printf("\t       B is P by N\n");
        printf("\t check!=0 executes also the standard ijk algo for checking the result\n\n");
        printf("\t NOTE: the blocksize must evenly divide M, N and P.\n");
        return -1;
    }

    long M = parse_arg(argv[1]);
    long N = parse_arg(argv[2]);
    long P = parse_arg(argv[3]);
    long b = parse_arg(argv[4]);
    long nw = parse_arg(argv[5]);
    if (argc >= 7) check = (atoi(argv[6])?true:false);

    const double *A = (double*)malloc(M*P*sizeof(double));
    const double *B = (double*)malloc(P*N*sizeof(double));
    assert(A); assert(B);

    random_init(M, N, P, const_cast<double*>(A), const_cast<double*>(B));
    double *C       = (double*)malloc(M*N*sizeof(double));

    ff_taskf taskf(nw);


    if (check) {
        double *C2   = (double*)malloc(M*N*sizeof(double));
        get_time(seqMatMult, NULL, 1, M,N,P, A, P, B, N, C, "Standard ijk algorithm");
        get_time_and_check(blockMMult, &taskf, b, M, N, P, A, P, B, N, C, C2, "Block-based MM");
        free(C2);
    } else 
        get_time(blockMMult, &taskf, b, M, N, P, A, P, B, N, C, "Block-based MM");
    
    free((void*)A); free((void*)B); free(C);
    
    return 0;  
} 
