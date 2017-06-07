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
 *  FastFlow version: 12/10/2009
 *
 *  August 2014: newer FastFlow version by Massimo Torquati
 *
 *  BUG TO BE FIXED: if thread workers deregister and re-register from the allocator several times, it crashes !!!
 *
 */

//#define ALLOCATOR_STATS 1 

#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <ff/allocator.hpp>
#include <ff/farm.hpp>
#include <ff/spin-lock.hpp>
#include "lib/marX2.h"

const int DIM=800;
const int ITERATION=1024;
const double init_a=-2.125, init_b=-1.5, range=3.0;   

using namespace ff;

static ff_allocator * ffalloc = NULL;
static lock_t lock;

struct ostream_t {
  unsigned char * M;
  int line;
};

struct Worker1: ff_node_t<ostream_t> {
    Worker1(const int dim, const int niter, const double step):
	dim(dim),niter(niter),step(step) {}

    ostream_t *svc(ostream_t *item) {
	int i = 	item->line; 
	double im=init_b+(step*i);
	for (int j=0;j<dim;j++)   {         
	    double a=init_a+step*j;
	    double b=im;
	    const double cr = a;
	    int k=0;
	    for ( ;k<niter;k++) {
		const double a2=a*a;
		const double b2=b*b;
		if ((a2+b2)>4.0) break;
		b=2*a*b+im;
		a=a2-b2+cr;
	    }
	    item->M[j] =  (unsigned char) 255-((k*255/niter)); 
	}  
	return item;	
    }

    const int dim,niter;
    const double step;
};

struct Worker2: ff_node_t<ostream_t> {
    Worker2(const int dim, const int niter, const double step):
	dim(dim),niter(niter),step(step) {}

    int svc_init() {
	if (ffalloc->register4free()<0) {
	    error("Worker2, register4free fails\n");
	    return -1;
	}
	return 0;
    }
    ostream_t *svc(ostream_t *item) {
	int i = 	item->line; 	  // this is just alpha renamining
	double im=init_b+(step*i);
	for (int j=0;j<dim;j++)   {         
	    double a=init_a+step*j;
	    double b=im;
	    const double cr = a;
	    int k=0;
	    for ( ;k<niter;k++) {
		const double a2=a*a;
		const double b2=b*b;
		if ((a2+b2)>4.0) break;
		b=2*a*b+im;
		a=a2-b2+cr;
	    }
	    item->M[j] =  (unsigned char) 255-((k*255/niter)); 
	  }        
#if !defined(NO_DISPLAY)
	spin_lock(lock);
	ShowLine(item->M,dim,item->line); 
	spin_unlock(lock);
#endif
	ffalloc->free(item->M);
	ffalloc->free(item);
	return GO_ON;
    }

    const int dim,niter;
    const double step;
};

struct Emitter: ff_node_t<ostream_t> {
public:
    Emitter(const int dim):
	dim(dim), ntask(dim) {}
    
    bool registered = false;

    int svc_init() {
	if (!registered)
	    if (ffalloc->registerAllocator()<0) {
		error("Emitter, registerAllocator fails\n");
		return -1;
	    }
	registered = true;
	ntask = dim;
	return 0;
    }
    ostream_t *svc(ostream_t *) {	
	int pos=dim-ntask;
	--ntask;
	if (ntask<0) return EOS;
	ostream_t * item = (ostream_t *)ffalloc->malloc(sizeof(ostream_t)); 
	item->M = (unsigned char *)ffalloc->malloc(dim*sizeof(char));
	item->line = pos;
	return item;
    }
    int dim,ntask;
};

struct Collector: ff_node_t<ostream_t> {
    Collector(const int dim):dim(dim) {}

    bool registered = false;

    int svc_init() {
	if (!registered) 
	    if (ffalloc->register4free()<0) {
		error("Worker, register4free fails\n");
		return -1;
	}
	return 0;
    }
    ostream_t *svc(ostream_t * t) {	  	
#if !defined(NO_DISPLAY)
	ShowLine(t->M,dim,t->line); 
#endif
	ffalloc->free(t->M);
	ffalloc->free(t);
	return GO_ON;
    }
    const int dim;
};


int main(int argc, char ** argv) {
    int ncores = ff_numCores();
    int dim = DIM, niter = ITERATION, retries=1, nworkers;

    if (argc<3) {
	printf("\nUsage: %s <nworkers> <size> [niterations=1024] [retries=1] [0|1]\n\n\n", argv[0]);
	return -1;
    } 
    nworkers = atoi(argv[1]);
    dim      = atoi(argv[2]);
    if (argc>=4) niter   = atoi(argv[3]);
    if (argc>=5) retries = atoi(argv[4]);
    // if argv[5] is defined and it is not zero, than we force to use the 
    // template with collector otherwise we need a lock to protect 
    // concurrent accesses to display lines
    if (argc>=6) { if (atoi(argv[5])) ncores=999; else ncores=2; }

    init_unlocked(lock);
    double step = range/((double)dim);
    double runs[retries];
    std::cout << "Mandebroot set from (" << init_a << "+I " << init_b << ") to (" 
	      << init_a+range << "+I " << init_b+range << ")\n";
    std::cout << "Resolution " << dim*dim << " pixels, max n. of iterations " << niter << "\n";
    
    if (ncores>=4) std::cout << "\nNOTE: using the farm template WITH the collector node!\n\n";
    else      	   std::cout << "\nNOTE: using the farm template WITHOUT the collector node!\n\n";
    
#if !defined(NO_DISPLAY)
    SetupXWindows(dim,dim,1,NULL,"FF Mandelbroot");
#endif

    std::vector<std::unique_ptr<ff_node> > w;
    for (int k=0;k<nworkers;k++)
	w.push_back((ncores>=4)? std::unique_ptr<ff_node>(new Worker1(dim,niter,step)) : 
		    std::unique_ptr<ff_node>(new Worker2(dim,niter,step)));
    Emitter E(dim);	
    Collector C(dim);
    ff_Farm<> farm(std::move(w),E,C); 
    if (ncores<4) {
	std::cout << "Collector removed!\n";
	farm.remove_collector();
    }

    // set auto-scheduling policy
    farm.set_scheduling_ondemand();

    // using the FastFlow allocator
    ffalloc = new ff_allocator;
    if (ffalloc && ffalloc->init()<0) {
	error("cannot build ff_allocator\n");
	return -1;
    }
        
    double avg=0.0, var=0.0;    
    for (int r=0;r<retries;r++) {
	ffTime(START_TIME);

	if (farm.run_then_freeze()<0) {
	    error("running farm\n");
	    return -1;
	}
	if (farm.wait_freezing()<0) {
	    error("wait freezing farm\n");
	    return -1;
	}

	ffTime(STOP_TIME);
	avg += runs[r] = ffTime(GET_TIME);
	
	std::cerr << "Run [" << r << "] DONE, time= " <<  runs[r] << " (ms)\n";	
    }
    delete ffalloc;
    avg /= (double) retries;
    for (int r=0;r<retries;r++) 
	var += (runs[r] - avg) * (runs[r] - avg);
    var /= (double) retries;
    std::cout << "The average time (ms) on " << retries << " experiments is " << avg 
	      << " Std. Dev. is " << sqrt(var) << "\n";
    
#if !defined(NO_DISPLAY)  
    std::cout << "Press a button to close the windows\n";    
    getchar();
    CloseXWindows();
#endif
    
    return 0;
}
