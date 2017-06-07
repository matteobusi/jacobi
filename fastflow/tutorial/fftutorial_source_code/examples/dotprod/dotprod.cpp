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
/*  This program computes the dot-product of 2 arrays as:
 *
 *    s = sum_i(A[i]*B[i])
 *
 *  command: 
 *
 *    dotprod ntimes size pfworkers G
 *
 *  ntimes the computation is repeated n times
 *  size is the size of the two arrays
 *  pfworkers is the number of workers threads used
 *  G=0 means default static scheduling
 *  G>0 means dynamic scheduling with grain N
 *  G<0 means static scheduling with grain N
 */

/* Author: Massimo Torquati
 * Date  : August 2014 
 *         
 */

#include <iostream>
#include <iomanip>
#include <ff/parallel_for.hpp>
using namespace ff;

int main(int argc, char * argv[]) {    
    if (argc < 5) {
        std::cerr << "use: " << argv[0] 
                  << " ntimes size pfworkers G\n";
        std::cerr << " example:\n  " << argv[0] << " 2 100000 8 0\n\n";
        return -1;
    }
    const double INITIAL_VALUE = 1.0;
    int  NTIMES    = atoi(argv[1]);
    long arraySize = atol(argv[2]);
    int  nworkers  = atoi(argv[3]);
    long chunk     = atol(argv[4]);

    assert(nworkers>0); assert(arraySize>0);
    
    // creates the array
    double *A = new double[arraySize];
    double *B = new double[arraySize];

    ParallelForReduce<double> pfr(nworkers,true); // spinwait is set to true
    auto Fsum = [](double& v, const double& elem) { v += elem; };
    
    double sum;
    if (chunk < 0) {  // static scheduling (round-robin with stride 'step')
        const long step = 1;
        pfr.parallel_for_static(0,arraySize, step, chunk, [&](const long j) { A[j]=j*3.14; B[j]=2.1*j;});

        ffTime(START_TIME);    
        for(int z=0;z<NTIMES;++z) {
            sum = INITIAL_VALUE;
            pfr.parallel_reduce_static(sum, 0.0,
                                       0, arraySize,1,chunk,
                                       [&](const long i, double& sum) { sum += A[i]*B[i]; }, 
                                       Fsum);
        }
        ffTime(STOP_TIME);
    } else {  /* dynamic scheduling with minumin granularity equal to chunk
               * note that if chunk=0 than it falls in the default static scheduling 
               * meaning that to each thread worker is assigned (statically) a big partition
               * of size ~(arraySize/nworkers).
               */
        pfr.parallel_for(0,arraySize,1, chunk, [&](const long j) { A[j]=j*3.14; B[j]=2.1*j;});

        ffTime(START_TIME);    
        for(int z=0;z<NTIMES;++z) {
            sum = INITIAL_VALUE;
            pfr.parallel_reduce(sum, 0.0, 
                                0, arraySize,1,chunk,
                                [&](const long i, double& sum) { sum += A[i]*B[i]; }, 
                                Fsum);
        }
        ffTime(STOP_TIME);
    }

    std::cout << "Result = " << std::setprecision(10) << sum << " Workers = " << nworkers << " Time (ms) = " << std::setprecision(5) << ffTime(GET_TIME)/(double)NTIMES << "\n";

    delete [] A;  delete [] B;
    return 0;
}
