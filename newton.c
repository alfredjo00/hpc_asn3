#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <threads.h>

typedef char TYPE_ATTR;
typedef char TYPE_CONV;

typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded;

typedef struct {
  const TYPE_ATTR **attractors;
  TYPE_CONV **convergences;
  int ib;
  int istep;
  int n_size;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  const TYPE_ATTR **attractors;
  TYPE_CONV **convergences;
  int n_size;
  int n_threads;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int
main_wrt_thrd(void *args){
    const thrd_info_t *thrd_info = (thrd_info_t*) args;
    const TYPE_ATTR **attractors = thrd_info->attractors;
    TYPE_CONV **convergences = thrd_info->convergences;
    const int istep = thrd_info->istep;
    const int n_size = thrd_info->n_size;
    const int ib = thrd_info->ib;

    printf("Writing thread: ib = %i\n", ib);

    return 0;
}

int
main_cmp_thrd(void *args){
    const thrd_info_t *thrd_info = (thrd_info_t*) args;
    const TYPE_ATTR **attractors = thrd_info->attractors;
    TYPE_CONV **convergences = thrd_info->convergences;
    const int istep = thrd_info->istep;
    const int n_size = thrd_info->n_size;
    const int ib = thrd_info->ib;

    TYPE_ATTR *attractor = (TYPE_ATTR*) malloc(n_size*n_size*sizeof(TYPE_ATTR));
    TYPE_CONV *convergence = (TYPE_CONV*) malloc(n_size*n_size*sizeof(TYPE_CONV));

    for ( size_t cx = 0; cx < n_size; ++cx ) {
        attractor[cx] = 0;
        convergence[cx] = 0;
    }

    printf("Computation thread: ib = %i\n", ib);
    
    return 0;
}

int
main_thrd(void* args){
    const thrd_info_t *thrd_info = (thrd_info_t*) args;
    const TYPE_ATTR **attractors = thrd_info->attractors;
    TYPE_CONV **convergences = thrd_info->convergences;
    const int istep = thrd_info->istep;
    const int n_size = thrd_info->n_size;
    const int ib = thrd_info->ib;

    printf("Main thread: ib = %i\n", ib);

    return 0;
}

int
main_thrd_check(void* args){
    const thrd_info_t *thrd_info = (thrd_info_t*) args;
    const TYPE_ATTR **attractors = thrd_info->attractors;
    TYPE_CONV **convergences = thrd_info->convergences;
    const int istep = thrd_info->istep;
    const int n_size = thrd_info->n_size;
    const int ib = thrd_info->ib;

    printf("Check thread: ib = %i\n", ib);
    
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
	
    TYPE_ATTR **attractors = (TYPE_ATTR**) malloc(n_size*sizeof(TYPE_ATTR*));
    TYPE_CONV **convergences = (TYPE_CONV**) malloc(n_size*sizeof(TYPE_CONV*));
    TYPE_ATTR *attractor = (TYPE_ATTR*) malloc(n_size*n_size*sizeof(TYPE_ATTR));

    thrd_t thrds[n_threads];
    thrd_info_t thrds_info[n_threads];

    thrd_t thrd_check;
    thrd_info_check_t thrd_info_check;
    
    mtx_t mtx;
    mtx_init(&mtx, mtx_plain);

    cnd_t cnd;
    cnd_init(&cnd);

    int_padded status[n_threads];

    // computation threads
    for ( int tx = 0; tx < n_threads; ++tx ) {
        thrds_info[tx].attractors = (const TYPE_ATTR**) attractors;
        thrds_info[tx].convergences = convergences;
        thrds_info[tx].ib = tx;
        thrds_info[tx].istep = n_threads;
        thrds_info[tx].n_size = n_size;
        thrds_info[tx].tx = tx;
        thrds_info[tx].mtx = &mtx;
        thrds_info[tx].cnd = &cnd;
        thrds_info[tx].status = status;
        status[tx].val = -1; // should this be zero?? check video

        int r = thrd_create(thrds+tx, main_cmp_thrd, (void*) (thrds_info+tx));
        if ( r != thrd_success ) {
        fprintf(stderr, "failed to create thread\n");
        exit(1);
        }
        thrd_detach(thrds[tx]);
    }

    // writing threads
    for ( int tx = 0; tx < n_threads; ++tx ) {
        thrds_info[tx].attractors = (const TYPE_ATTR**) attractors;
        thrds_info[tx].convergences = convergences;
        thrds_info[tx].ib = tx;
        thrds_info[tx].istep = n_threads;
        thrds_info[tx].n_size = n_size;
        thrds_info[tx].tx = tx;
        thrds_info[tx].mtx = &mtx;
        thrds_info[tx].cnd = &cnd;
        thrds_info[tx].status = status;
        status[tx].val = -1; // should this be zero?? check video

        int r = thrd_create(thrds+tx, main_wrt_thrd, (void*) (thrds_info+tx));
        if ( r != thrd_success ) {
        fprintf(stderr, "failed to create thread\n");
        exit(1);
        }
        thrd_detach(thrds[tx]);
    }

    {
        thrd_info_check.attractors = (const TYPE_ATTR**) attractors;
        thrd_info_check.convergences = convergences;
        thrd_info_check.n_size = n_size;
        thrd_info_check.n_threads = n_threads;
        thrd_info_check.mtx = &mtx;
        thrd_info_check.cnd = &cnd;
        // It is important that we have initialize status in the previous for-loop,
        // since it will be consumed by the check threads.
        thrd_info_check.status = status;

        int r = thrd_create(&thrd_check, main_thrd_check, (void*) (&thrd_info_check));
        if ( r != thrd_success ) {
        fprintf(stderr, "failed to create thread\n");
        exit(1);
        }
    }

    {
        int r;
        thrd_join(thrd_check, &r);
    }

    free(attractor);
    free(attractors);
    free(convergences);
	
   return 0;
}