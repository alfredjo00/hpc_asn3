#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <threads.h>

typedef struct {
  const float **v;
  float **w;
  int ib;
  int istep;
  int sz;
  int tx;
  mtx_t *mtx;
} thrd_info_t;

int thrd_fun(void *args)
{
	const thrd_info_t *thrd_info = (thrd_info_t*) args;
	const float **v = thrd_info->v;
	float **w = thrd_info->w;
	const int ib = thrd_info->ib;
	const int istep = thrd_info->istep;
	const int sz = thrd_info->sz;
	const int tx = thrd_info->tx;
	mtx_t *mtx = thrd_info->mtx;
	
	for ( int ix = ib; ix < sz; ix += istep ) {
		const float *vix = v[ix];
		// We allocate the rows of the result before computing, and free them in another thread.
		float *wix = (float*) malloc(3*sz*sizeof(float));

		for ( int jx = 0; jx < sz; ++jx ){			
			wix[3 * jx] 		= jx % 20; 	// R
			wix[3 * jx + 1] 	= jx % 100; // G
			wix[3 * jx + 2] 	= jx % 200; // B
		}

		mtx_lock(mtx);
		w[ix] = wix;
		mtx_unlock(mtx);
	}
	printf("Done %d\n", tx);
	return 0;
}

int write_ppm(float **w, int n_size)
{
	
	
	for (int i = 0; i < n_size; i++)
		free(w[i]);
}


int main(int argc, char *argv[]) {
   int n_size = 0, n_threads = 0;

   for (int i = 1; i < argc; i++) {		
      if (strncmp(argv[i], "-l", 2) == 0) 
         n_size = atoi(argv[i] + 2);
		
      else if (strncmp(argv[i], "-t", 2) == 0) 
         n_threads = atoi(argv[i] + 2);      
   }
	
	assert(n_size != 0 && n_threads != 0);	
	
	float **w = (float**) malloc(n_size*sizeof(float*));
	float **v = (float**) malloc(n_size*sizeof(float*));
	float *ventries = (float*) malloc(2*n_size*n_size*sizeof(float));


	for ( int ix = 0, jx = 0; ix < n_size; ++ix, jx += 2 * n_size )
		v[ix] = ventries + jx;

	float delta = 4 / n_size;
	for ( int ix = 0; ix < n_size; ++ix )
	{
		for ( int jx = 0; jx < n_size; ++jx )
		{	
			// Real and imaginary coords. (re, im)
			ventries[2 * ix * n_size + 2 * jx] 		= -2 + ix * delta;
			ventries[2 * ix * n_size + 2 * jx + 1] =  2 - jx * delta;
		}
	}
	
	thrd_t thrds[n_threads];
	thrd_info_t thrds_info[n_threads];

	mtx_t mtx;
	mtx_init(&mtx, mtx_plain);

	for ( int tx = 0; tx < n_threads; ++tx ) {
		thrds_info[tx].v = (const float**) v;
		thrds_info[tx].w = w;
		thrds_info[tx].ib = tx;
		thrds_info[tx].istep = n_threads;
		thrds_info[tx].sz = n_size;
		thrds_info[tx].tx = tx;
		thrds_info[tx].mtx = &mtx;

		int r = thrd_create(thrds+tx, thrd_fun, (void*) (thrds_info+tx));
		if ( r != thrd_success ) {
			fprintf(stderr, "failed to create thread\n");
			exit(1);
		}
	}	

	for ( int tx = 0; tx < n_threads; ++tx ) {
		int r;
		thrd_join(thrds[tx], &r);
	}

	// Write to ppm files

	mtx_destroy(&mtx);
   return 0;
}
