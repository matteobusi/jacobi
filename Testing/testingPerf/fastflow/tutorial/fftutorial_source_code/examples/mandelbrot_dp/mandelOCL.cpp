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
 *  August 2014: Massimo Torquati
 *                    
 */

#if !defined(FF_OPENCL)
#define FF_OPENCL
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <ff/utils.hpp>
#include <ff/stencilReduceOCL.hpp>
#include "lib/marX2.h"    

const int DIM=800;
const int ITERATION=1024;
const double init_a=-2.125, init_b=-1.5, range=3.0;   

using namespace ff;

struct parameters_t {
    int   niter;
    float init_a;
    float init_b;
    float step;
};
#if !defined(BUILD_WITH_SOURCE)
FF_OCL_STENCIL_ELEMFUNC_2D_IO(mapf, parameters_t, unsigned char, h,w,i,j,
			      float  im=in->init_b+(in->step*i);
			      float  a =in->init_a+(in->step*j);
			      float  b=im;
			      const float cr = a;
			      int k=0;
			      //printf("%d, %d, im=%f, a=%f, b=%f, niter=%d\n", i,j,im, a, b, in->niter);
			      for ( ; k<in->niter;k++) {
				  const float  a2=a*a;
				  const float b2=b*b;
				  if ((a2+b2)>4.0) break;
				  b=2*a*b+im;
				  a=a2-b2+cr;
			      }
			      //printf("k=%d, %d\n", k, (unsigned char) (255-((k*255/in->niter))) );
			      return (unsigned char) (255-((k*255/in->niter))); 
);
#endif

struct oclTaskM: public baseOCLTask_2D<oclTaskM, parameters_t, unsigned char> {
    oclTaskM():p(nullptr),M(nullptr),dimX(0),dimY(0) {}

    oclTaskM(parameters_t *p, unsigned char *M, int dimX, int dimY):
	p(p),M(M),dimX(dimX),dimY(dimY) {}

    void setTask(const oclTaskM *t) {
	setHeight(t->dimX);
	setWidth(t->dimY);
	setInPtr(t->p, 1);
	setOutPtr(t->M, t->dimX*t->dimY);
    }

    parameters_t *p;
    unsigned char *M;
    int dimX;
    int dimY;
};



int main(int argc, char **argv) {
    int dim = DIM, niter = ITERATION;

    if (argc >=2) dim   = atoi(argv[1]);
    if (argc >=3) niter = atoi(argv[2]);

    double step = range/((double) dim);
    
    std::cout << "Mandebroot set from (" << init_a << "+I " << init_b << ") to (" 
	      << init_a+range << "+I " << init_b+range << ")\n";
    std::cout << "Resolution " << dim*dim << " pixels, max n. of iterations " << niter << "\n";
    
#if !defined(NO_DISPLAY)
    SetupXWindows(dim,dim,1,NULL,"ParallelForPipeReduce Mandelbrot");
#endif

    int device;
    std::vector<cl_device_id> dev;
    if ((device = clEnvironment::instance()->getGPUDevice()) == -1) {
	std::cerr << "NO GPU device available, switching to CPU\n";
	device = clEnvironment::instance()->getCPUDevice(); // this exists for sure
    } 
    dev.push_back(clEnvironment::instance()->getDevice(device)); //convert logical to OpenCL Id
        
    parameters_t param;
    param.niter  = niter;
    param.init_a = init_a;
    param.init_b = init_b;
    param.step   = step;

    std::vector<unsigned char> M(dim*dim, 0);
    oclTaskM oclt(&param, const_cast<unsigned char*>(M.data()), dim, dim);
    
#if defined(BUILD_WITH_SOURCE)
    ff_mapOCL_2D<oclTaskM> map2D(oclt, "code_cl/mandel.cl", "mapf");
    map2D.saveBinaryFile(); 
    map2D.reuseBinaryFile();
#else
    ff_mapOCL_2D<oclTaskM> map2D(oclt, mapf);
#endif
    map2D.setDevices(dev);
    
    ffTime(START_TIME);
    if (map2D.run_and_wait_end()<0) {
	error("running map\n");
	return -1;
    }
    ffTime(STOP_TIME);
    printf("Time= %f (ms)\n", ffTime(GET_TIME));
        
#if !defined(NO_DISPLAY)
    for(int i=0;i<dim;++i)
	ShowLine(&M[i*dim],dim,i); 

    std::cout << "Press a button to close the windows\n";
    getchar();
    CloseXWindows();
#endif   

    return 0;
}
