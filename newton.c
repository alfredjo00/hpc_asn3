#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <threads.h>
#include <math.h>
#include <errno.h>
#include <limits.h>
	

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))
#define QD(x) ((x)*(x)*(x)*(x))
#define QN(x) ((x)*(x)*(x)*(x)*(x))


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
	int d;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  float **w;
	int **f;
  int sz;
	int d;
  int nthrds;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;


static void poly_compute(float x, float y, float* z, int d)
{
	switch(d)
	{
		case 1:
		{
			z[0] = x - 1.0f;
			z[1] = y - 1.0f;
			break;
		}
		case 2:
		{
			float sq_x = SQ(x), sq_y = SQ(y);
			z[0] = x * (sq_x + sq_y - 1.0f) / (2.0f * ( sq_x + sq_y));
			z[1] = y * (sq_x + sq_y + 1.0f) / (2.0f * ( sq_x + sq_y));
			break;
		}
		case 3:
		{
			float sq_x = SQ(x), sq_y = SQ(y);
			float cu_x = CU(x), cu_y = CU(y);
			float qd_x = QD(x), qd_y = QD(y);
			float qn_x = QN(x);
			
			float c = 3.0f * SQ(sq_x + sq_y);
			z[0] = (qn_x + 2.0f * cu_x * sq_y - sq_x + x * qd_y + sq_y) / c;
			z[1] = y * (qd_x + 2.0f * sq_x * sq_y + 2.0f*x + qd_y) / c;
			break;
		}
		case 4:
		{
			float sq_x = SQ(x), sq_y = SQ(y);
			float cu_x = CU(x), cu_y = CU(y);
			float qd_x = QD(x), qd_y = QD(y);
			float qn_x = QN(x), qn_y = QN(y);
			float he_x = SQ(cu_x), he_y = SQ(cu_y);
			float sv_x = x*he_x, sv_y = y*he_y;
			float c = 4.0f * CU(sq_x + sq_y);
			z[0] =  sv_x + 3.0f*qn_x*sq_y + cu_x*(3.0f*qd_y - 1.0f) + x*sq_y*(qd_y + 3.0f);
			z[1] = (he_x + 3.0f*qd_x*sq_y + 3.0f*sq_x*(qd_y + 1.0f) + sq_y*(qd_y - 1.0f)) * y / c;
			break;
		}
		case 5:
		{
			float sq_x = SQ(x), sq_y = SQ(y);
			float cu_x = CU(x), cu_y = CU(y);
			float qd_x = QD(x), qd_y = QD(y);
			float qn_x = QN(x), qn_y = QN(y);
			float he_x = SQ(cu_x), he_y = SQ(cu_y);
			float sv_x = x*he_x, sv_y = y*he_y;
			
			float c = 5.0f * QD(sq_x + sq_y);
			z[0] = (sq_x*sv_x + 4.0f*sv_x*sq_y + 6.0f*qn_x*qd_y - qd_x + 4.0f*cu_x*he_y +\
							6.0f*sq_x*sq_y + x*y*sv_y - qd_y)/c;
							
			z[1] = (x*sv_x + 4.0f*he_x*sq_y + 6.0f*qd_x*qd_y + 4.0f*cu_x + 4.0f*qd_x*he_y -\
							4.0f*x*sq_y + y*sv_y) * y/c;
			break;
		}
			
		default:
			fprintf(stderr, "unexpected degree\n");
			exit(1);
	}	
}

static int poly_iteration(float x, float y, int d)
{
	float roots_x[d];
	float roots_y[d];
	float pi = 3.141596536;
	for (int r = 0; r < d; r++)
	{
		roots_x[r] = cos(r * 2.0f * pi/d);
		roots_y[r] = sin(r * 2.0f * pi/d);	
	}
	float z[2];
	float zr = x , zi = y;
	float sq_norm, root_sq_norm;
	for (int i = 0; i < 100; i++){
		poly_compute(zr, zi, z, d);
		zr += -z[0];
		zi += -z[1];
		
		if (zr > 1e10 || zi > 1e10)
			return i;
		sq_norm = x*x + y*y;
		
		if (sq_norm < 1e-6)
			return i;
		
		for (int j = 0; j < d; j++){
			root_sq_norm = SQ(zr - roots_x[j]) + SQ(zi - roots_y[j]);
			if (root_sq_norm < 1e-6)
				return i;
		}
	}	
}

int thrd_fun(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  float **w 					= thrd_info->w;
  int **f 						= thrd_info->f;
  const int ib 				= thrd_info->ib;
  const int istep			= thrd_info->istep;
  const int sz 				= thrd_info->sz;
  const int d 				= thrd_info->d;
  const int tx 				= thrd_info->tx;
  mtx_t *mtx 					= thrd_info->mtx;
  cnd_t *cnd 					= thrd_info->cnd;
  int_padded *status 	= thrd_info->status;
	const float dxy   	= 4.0 /sz;
	short iters;
	for ( int ix = ib; ix < sz; ix += istep ) {
		// We allocate the rows of the result before computing, and free them in another thread.
		float *wix 	= (float*) malloc(3*sz*sizeof(float));
		int *fix 		= (int*) 	malloc(sz*sizeof(int));
		float z[2];
		
		for ( int jx = 0; jx < sz; ++jx ){		
			iters = poly_iteration(-2.0f + dxy*ix, 2.0f - dxy*jx, d);
			wix[3 * jx] 			= jx % 20; 	// R
			wix[3 * jx + 1] 	= jx % 100; // G
			wix[3 * jx + 2] 	= jx % 200; // B
			fix[jx] 					= iters;
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
  const int d 				= thrd_info->d;
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
	int n_size = 0, n_threads = 0, n_d = 0;

	for (int i = 1; i < argc; i++) {		
		if (strncmp(argv[i], "-l", 2) == 0) 
			n_size = atoi(argv[i] + 2);

		else if (strncmp(argv[i], "-t", 2) == 0) 
			n_threads = atoi(argv[i] + 2);  
	}
	 
	{
		char *p;
		errno = 0;
		n_d = strtol(argv[argc - 1], &p, 10);		
		
		if (errno != 0 || *p != '\0' || n_d > INT_MAX || n_d < INT_MIN)
			{
				printf("error with args.\n");
				exit(1);
			}
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
    thrds_info[tx].d 			= n_d;
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
    thrd_info_check.d 			= n_d;
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
