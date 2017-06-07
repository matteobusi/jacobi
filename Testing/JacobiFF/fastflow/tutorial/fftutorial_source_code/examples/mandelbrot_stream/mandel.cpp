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
 *  Author: Marco Aldinucci <aldinuc@di.unito.it>
 *  Date :  15/11/97
 *         
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <ff/utils.hpp>
#include "lib/marX2.h"    

const int DIM=800;
const int ITERATION=1024;
const double init_a=-2.125, init_b=-1.5, range=3.0;   

using namespace ff;

int main(int argc, char **argv) {
    int dim = DIM, niter = ITERATION, retries=1;
        
    if (argc<2) {
	printf("Usage: %s <size> [niterations] [retries]\n\n\n", argv[0]);
	return -1;
    } 
    dim = atoi(argv[1]);
    if (argc>=3) niter   = atoi(argv[2]);
    if (argc>=4) retries = atoi(argv[3]);

    double step = range/((double) dim);
    double runs[retries];
    unsigned char *M = (unsigned char *) malloc(dim*sizeof(char));
    assert(M);
    
    std::cout << "Mandebroot set from (" << init_a << "+I " << init_b << ") to (" 
	      << init_a+range << "+I " << init_b+range << ")\n";
    std::cout << "Resolution " << dim*dim << " pixels, max n. of iterations " << niter << "\n";
    
#if !defined(NO_DISPLAY)
    SetupXWindows(dim,dim,1,NULL,"Sequential Mandelbroot");
#endif

    double avg=0.0, var=0.0;    
    for (int r=0;r<retries;r++) {
	ffTime(START_TIME);
	
	for(int i=0;i<dim;i++) {
	    double im=init_b+(step*i);
	    for (int j=0;j<dim;j++) {         
		double a=init_a+step*j;
		double b=im;
		const double cr = a;
		int k=0;
		for ( ; k<niter;k++) {
		    const double a2=a*a;
		    const double b2=b*b;
		    if ((a2+b2)>4.0) break;
		    b=2*a*b+im;
		    a=a2-b2+cr;
		}
		M[j]= (unsigned char) 255-((k*255/niter)); 
	    }
#if !defined(NO_DISPLAY)
	    ShowLine(M,dim,i); 
#endif
	}
	ffTime(STOP_TIME);
	avg += runs[r] = ffTime(GET_TIME);
	std::cout << "Run [" << r << "] done, time = " << runs[r] << "\n";
    }
    avg = avg / (double) retries;
    for (int r=0;r<retries;r++) 
	var += (runs[r] - avg) * (runs[r] - avg);
    var /= retries;
    std::cout << "The average time (ms) on " << retries << " experiments is " << avg 
	      << " Std. Dev. is " << sqrt(var) << "\n";
        
#if !defined(NO_DISPLAY)
    std::cout << "Press a button to close the windows\n";
    getchar();
    CloseXWindows();
#endif   

    return 0;
}
