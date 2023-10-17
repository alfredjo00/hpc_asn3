#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <threads.h>


typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded;

typedef struct {
  float **w;
	int **f;
  int ib;
  int istep;
  int sz;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  float **w;
	int **f;
  int sz;
  int nthrds;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;


int thrd_fun(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  float **w 					= thrd_info->w;
  int **f 						= thrd_info->f;
  const int ib 				= thrd_info->ib;
  const int istep			= thrd_info->istep;
  const int sz 				= thrd_info->sz;
  const int tx 				= thrd_info->tx;
  mtx_t *mtx 					= thrd_info->mtx;
  cnd_t *cnd 					= thrd_info->cnd;
  int_padded *status 	= thrd_info->status;
	
	for ( int ix = ib; ix < sz; ix += istep ) {
		// We allocate the rows of the result before computing, and free them in another thread.
		float *wix 	= (float*) malloc(3*sz*sizeof(float));
		int *fix 		= (int*) 	malloc(sz*sizeof(int));

		for ( int jx = 0; jx < sz; ++jx ){			
			wix[3 * jx] 			= jx % 20; 	// R
			wix[3 * jx + 1] 	= jx % 100; // G
			wix[3 * jx + 2] 	= jx % 200; // B
			fix[jx] 					= jx % 20;
		}

		mtx_lock(mtx);
		w[ix] = wix; // array of RGB
		f[ix] = fix; // array of iterations done
    status[tx].val = ix + istep;
    mtx_unlock(mtx);
    cnd_signal(cnd);
	}
	printf("Done %d\n", tx);
	return 0;
}


int thrd_check_fun(void *args)
{
  const thrd_info_check_t *thrd_info = (thrd_info_check_t*) args;
  float **w 					= thrd_info->w;
  int **f 						= thrd_info->f;
  const int sz 				= thrd_info->sz;
  const int nthrds 		= thrd_info->nthrds;
  mtx_t *mtx 					= thrd_info->mtx;
  cnd_t *cnd 					= thrd_info->cnd;
  int_padded *status 	= thrd_info->status;

  // We do not increment ix in this loop, but in the inner one.
  for ( int ix = 0, ibnd; ix < sz; ) {

    // If no new lines are available, we wait.
    for ( mtx_lock(mtx); ; ) {
      // We extract the minimum of all status variables.
      ibnd = sz;
      for ( int tx = 0; tx < nthrds; ++tx )
        if ( ibnd > status[tx].val )
          ibnd = status[tx].val;

      if ( ibnd <= ix )
        // We rely on spurious wake-ups, which in practice happen, but are not
        // guaranteed.
        cnd_wait(cnd,mtx);
      else {
        mtx_unlock(mtx);
        break;
      }
    }

    // We do not initialize ix in this loop, but in the outer one.
    for ( ; ix < ibnd; ++ix ) {
      // We free the component of w, since it will never be used again.
      free(w[ix]);
      free(f[ix]);
    }
  }

  return 0;
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
	int **f 	= (int**) malloc(n_size*sizeof(int*));

  thrd_t thrds[n_threads];
  thrd_info_t thrds_info[n_threads];

  thrd_t thrd_check;
  thrd_info_check_t thrd_info_check;
  
  mtx_t mtx;
  mtx_init(&mtx, mtx_plain);

  cnd_t cnd;
  cnd_init(&cnd);

  int_padded status[n_threads];

	for ( int tx = 0; tx < n_threads; ++tx ) {
    thrds_info[tx].w 			= w;
    thrds_info[tx].f 			= f;
    thrds_info[tx].ib	 		= tx;
    thrds_info[tx].istep 	= n_threads;
    thrds_info[tx].sz 		= n_size;
    thrds_info[tx].tx 		= tx;
    thrds_info[tx].mtx 		= &mtx;
    thrds_info[tx].cnd 		= &cnd;
    thrds_info[tx].status = status;
    status[tx].val 				= 0;

		int r = thrd_create(thrds+tx, thrd_fun, (void*) (thrds_info+tx));
		if ( r != thrd_success ) {
			fprintf(stderr, "failed to create thread\n");
			exit(1);
		}
    thrd_detach(thrds[tx]);
	}	

  {
    thrd_info_check.w 			= w;
    thrd_info_check.f 			= f;
    thrd_info_check.sz 			= n_size;
    thrd_info_check.nthrds 	= n_threads;
    thrd_info_check.mtx 		= &mtx;
    thrd_info_check.cnd 		= &cnd;
    // It is important that we have initialize status in the previous for-loop,
    // since it will be consumed by the check threads.
    thrd_info_check.status = status;

    int r = thrd_create(&thrd_check, thrd_check_fun, (void*) (&thrd_info_check));
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }

  {
    int r;
    thrd_join(thrd_check, &r);
  }

	// Write to ppm files

  free(w);

  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
	
	return 0;
}
