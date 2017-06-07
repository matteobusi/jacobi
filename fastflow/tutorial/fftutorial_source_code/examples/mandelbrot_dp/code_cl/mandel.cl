typedef struct _parameters_t {
    int   niter;
    float init_a;
    float init_b;
    float step;
} parameters_t;


#define GET_IN(i,j) (in[((i)*w+(j))-offset])
#define GET_ENV1(i,j) (env1[((i)*w+(j))])
unsigned char fmapf(
	__global parameters_t* in,
	const uint h,
	const uint w,
	const int i,
	const int j,
	const int offset) {
	   float im=in->init_b+(in->step*i); 
	   float a =in->init_a+(in->step*j); 
	   float b=im; 
	   const float cr = a; 
	   int k=0; 
	   for ( ; k<in->niter;k++) { 
	     const float a2=a*a; 
	     const float b2=b*b; 
	     if ((a2+b2)>4.0) break; 
	     b=2*a*b+im; 
	     a=a2-b2+cr; 
	   } 
	   return (unsigned char) (255-((k*255/in->niter)));
}

__kernel void kern_mapf(
	__global parameters_t* input,
	__global unsigned char* output,
	const uint inHeight,
	const uint inWidth,
	const uint maxItems,
	const uint offset,
	const uint halo) {
	    size_t i = get_global_id(0);
	    size_t ig = i + offset;
	    size_t r = ig / inWidth;
	    size_t c = ig % inWidth;
	    size_t gridSize = get_local_size(0)*get_num_groups(0);
	    while(i < maxItems)  {
	        output[i+halo] = fmapf(input+halo,inHeight,inWidth,r,c,offset);
	        i += gridSize;
	    }
}
