#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <threads.h>
#include <math.h>
#include <errno.h>
#include <limits.h>	

#define SQ(x) ((x)*(x))

#define MIN(a,b) ({ __typeof__ (a) _a = (a); \
					__typeof__ (b) _b = (b); \
								_a < _b ? _a : _b; })

typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded;

typedef struct {
	int **attractors;
	int **convergences;
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
  	int **attractors;
	int **convergences;
  	int max_col_val;
  	FILE *file_conv;
  	FILE *file_attr;
  	int sz;
	int d;
  	int nthrds;
  	mtx_t *mtx;
  	cnd_t *cnd;
  	int_padded *status;
} thrd_info_check_t;

static void
write_header(FILE *file, int n_size, int max_color_val){
    fprintf(file, "%s\n%i %i\n%i\n", "P3", n_size, n_size, max_color_val);
}

static void
write_conv(FILE *file, int* convergence, int n_size, char colors[], char *row_str)
{
  int color_str_len = 20;
  int row_str_len_sum = 0;
  int offset = 0;

  // memcpy color strings to row_str
  for ( int ix = 0, jx = 0; jx < n_size; ix += color_str_len, ++jx ){
    assert( 1 <= convergence[jx] && convergence[jx] <= 128 );
      if ( convergence[jx] <= 5 ){
          color_str_len = 6;
          row_str_len_sum += color_str_len;
          offset = 0;
      }
      else if ( convergence[jx] <= 50 ){
          color_str_len = 9;
          row_str_len_sum += color_str_len;
          offset = -3*5;
      }
      else {
          color_str_len = 12;
          row_str_len_sum += color_str_len;
          offset = -6*5 - 3*45;
      }
      memcpy( row_str + ix, colors + offset + color_str_len*(convergence[jx]-1), color_str_len);
  }
  fwrite(row_str, sizeof(char), row_str_len_sum, file);
}

static void
write_attr(FILE *file, int* attractor, int n_size, int n_degree, char colors[], char* row_str)
{
    int color_str_len = 12;

    for ( size_t ix = 0, jx = 0; jx < n_size; ix += color_str_len, ++jx ){
      	assert( attractor[jx] <= 10 && attractor[jx] >= 0);
      	memcpy( row_str + ix, colors + attractor[jx]*color_str_len, color_str_len);        
    }
    fwrite( row_str, sizeof(char), n_size*color_str_len, file);
}

static void poly_compute(float x, float y, float* z, int d)
{
	switch(d)
	{
		case 1:
		{
			z[0] = x - 1.0f;
			z[1] = y;
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
			float sq_x = SQ(x), 	sq_y = SQ(y);
			float cu_x = x*sq_x, 	cu_y = y*sq_y;
			float qd_x = x*cu_x, 	qd_y = y*cu_y;
			float qn_x = x*qd_x;
			
			float c = 3.0f * SQ(sq_x + sq_y);
			z[0] = (qn_x + 2.0f * cu_x * sq_y - sq_x + x * qd_y + sq_y) / c;
			z[1] = y * (qd_x + 2.0f * sq_x * sq_y + 2.0f*x + qd_y) / c;
			break;
		}
		case 4:
		{
			float sq_x = SQ(x), 	sq_y = SQ(y);
			float cu_x = x*sq_x, 	cu_y = y*sq_y;
			float qd_x = x*cu_x, 	qd_y = y*cu_y;
			float qn_x = x*qd_x,	qn_y = y*qd_y;
			float he_x = x*qn_x,	he_y = y*qn_y;
			float sv_x = x*he_x,	sv_y = y*he_y;
			float oc_x = x*sv_x,	oc_y = y*sv_y;
			
			float c = 4.0f*he_x + 12.0f*qd_x*sq_y + 12.0f*sq_x*qd_y + 4.0f*he_y;
			z[0] = sv_x + 3.0f*qn_x*sq_y + cu_x*(3.0f*qd_y - 1.0f) + x*sq_y*(qd_y + 3.0f);
			z[0] = z[0]/c;
			
			z[1] = (he_x + 3.0f*qd_x*sq_y + 3.0f*sq_x*(qd_y + 1.0f) + sq_y*(qd_y - 1.0f)) * y;
			z[1] = z[1]/c;
			break;
		}
		case 5:
		{
			float sq_x = SQ(x), 	sq_y = SQ(y);
			float cu_x = x*sq_x, 	cu_y = y*sq_y;
			float qd_x = x*cu_x, 	qd_y = y*cu_y;
			float qn_x = x*qd_x,	qn_y = y*qd_y;
			float he_x = x*qn_x,	he_y = y*qn_y;
			float sv_x = x*he_x,	sv_y = y*he_y;
			float oc_x = x*sv_x,	oc_y = y*sv_y;
			
			
			float c = 5.0f*oc_x + 20.0f*he_x*sq_y + 30.0f*qd_x*qd_y + 20.0f*sq_x*he_y + 5.0f*oc_y;
			float a = oc_x + 4.0f*he_x*sq_y + 6.0f*qd_x*qd_y + 4.0f*sq_x*he_y + oc_y;
			z[0] = x*a - qd_x + 6.0f*sq_x*sq_y - qd_y;
			z[0] = z[0]/c;

			z[1] = y*(a + 4.0f*cu_x - 4.0f*x*sq_y);
			z[1] = z[1]/c;
			break;
		}
			
		case 6:
		{
			float sq_x = SQ(x), 		sq_y = SQ(y);
			float cu_x = x*sq_x, 		cu_y = y*sq_y;
			float qd_x = x*cu_x, 		qd_y = y*cu_y;
			float qn_x = x*qd_x,		qn_y = y*qd_y;
			float he_x = x*qn_x,		he_y = y*qn_y;
			float sv_x = x*he_x,		sv_y = y*he_y;
			float oc_x = x*sv_x,		oc_y = y*sv_y;
			float tn_x = sq_x*oc_x, tn_y = sq_y*oc_y;
			
			
			float c = 6.0f*tn_x + 6.0f*tn_y + 30.0f*oc_x*sq_y + 30.0f*oc_y*sq_x + 60.0f*he_x*qd_y + 60.0f*he_y*qd_x;
			float a = tn_x + 5.*oc_x*sq_y + 10.*he_x*qd_y + 5.*sq_x*oc_y;
			
			z[0] = x*a + qn_x*(10.0f*he_y-1.0f) + 10.0f*cu_x*sq_y + x*tn_y - 5.0f*x*qd_y;
			z[0] = z[0]/c;

			z[1] = y*(a + 5.0f*qd_x + 10.0f*qd_x*he_y - 10*sq_x*sq_y + tn_y + qd_y);
			z[1] = z[1]/c;
			break;
		}
			
		case 7:
		{
			float sq_x = SQ(x), 		sq_y = SQ(y);
			float cu_x = x*sq_x, 		cu_y = y*sq_y;
			float qd_x = x*cu_x, 		qd_y = y*cu_y;
			float qn_x = x*qd_x,		qn_y = y*qd_y;
			float he_x = x*qn_x,		he_y = y*qn_y;
			float sv_x = x*he_x,		sv_y = y*he_y;
			float oc_x = x*sv_x,		oc_y = y*sv_y;
			float tn_x = sq_x*oc_x, tn_y = sq_y*oc_y;
			float tw_x = sq_x*tn_x, tw_y = sq_y*tn_y;
			
			
			float c = 7.0f*(tw_x + tw_y) + 42.0f*(tn_x*sq_y + tn_y*sq_x) +\
								105.0f*(oc_x*qd_y + oc_y*qd_x) + 140.0f*he_x*he_y;
								
			float a = tw_x + 6.0f*tn_x*sq_y + 6.0f*sq_x*tn_y + 15.0f*oc_x*qd_y + 15.0f*qd_x*oc_y +\
								tw_y + 20.0f*he_x*he_y;
			
			z[0] = x*a - he_x + 15.0f*(qd_x*sq_y - sq_x*qd_y) + he_y;
			z[0] = z[0]/c;
				
			z[1] = y*a + 6.0f*(x*qn_y + y*qn_x) - 20.0f*cu_x*cu_y;
			z[1] = z[1]/c;
			break;
		}		
				
		case 8:
		{
			float sq_x = SQ(x), 		sq_y = SQ(y);
			float cu_x = x*sq_x, 		cu_y = y*sq_y;
			float qd_x = x*cu_x, 		qd_y = y*cu_y;
			float qn_x = x*qd_x,		qn_y = y*qd_y;
			float he_x = x*qn_x,		he_y = y*qn_y;
			float sv_x = x*he_x,		sv_y = y*he_y;
			float oc_x = x*sv_x,		oc_y = y*sv_y;
			float tn_x = sq_x*oc_x, tn_y = sq_y*oc_y;
			float tw_x = sq_x*tn_x, tw_y = sq_y*tn_y;
			float ft_x = sq_x*tw_x, ft_y = sq_y*tw_y;
			
			float c = 8.0f*(ft_x + ft_y) + 56.0f*(tw_x*sq_y + tw_y*sq_x) +\
								168.0f*(tn_x*qd_y + tn_y*qd_x) + 280.0f*(oc_x*he_y + oc_y*he_x);
								
			float a = ft_y + ft_x + 7.0f*(sq_x*tw_y + tw_x*sq_y) + 21.0f*(qd_x*tn_y + tn_x*qd_y) +\
								35.0f*(oc_x*he_y + oc_y*he_x);			
						
			z[0] = x*a - sv_x - 35.0f*cu_x*qd_y + 21.0f*qn_x*sq_y + 7.0f*x*he_y;
			z[0] = z[0]/c;

			z[1] = y*a - sv_y - 35.0f*qd_x*cu_y + 21.0f*sq_x*qn_y + 7.0f*he_x*y;
			z[1] = z[1]/c;
			break;
		}				

		case 9:
		{
			float sq_x = SQ(x), 		sq_y = SQ(y);
			float cu_x = x*sq_x, 		cu_y = y*sq_y;
			float qd_x = x*cu_x, 		qd_y = y*cu_y;
			float qn_x = x*qd_x,		qn_y = y*qd_y;
			float he_x = x*qn_x,		he_y = y*qn_y;
			float sv_x = x*he_x,		sv_y = y*he_y;
			float oc_x = x*sv_x,		oc_y = y*sv_y;
			float tn_x = sq_x*oc_x, tn_y = sq_y*oc_y;
			float tw_x = sq_x*tn_x, tw_y = sq_y*tn_y;
			float ft_x = sq_x*tw_x, ft_y = sq_y*tw_y;
			float sx_x = sq_x*ft_x, sx_y = sq_y*ft_y;			
			
			float c = 9.0f*(sx_x + sx_y) + 72.0f*(ft_x*sq_y + ft_y*sq_x) + 252.0f*(tw_x*qd_y + tw_y*qd_x) + 504.0f*(tn_x*he_y + tn_y*he_x) + 630.0f*oc_x*oc_y;
			float a = sx_x + sx_y + 8.0f*(ft_x*sq_y + ft_y*sq_x) + 28.0f*(tw_x*qd_y + tw_y*qd_x) + 56.0f*(tn_x*he_y + tn_y*he_x) + 70.0f*oc_x*oc_y;
			
			z[0] = x*a - oc_x - oc_y + 28.0f*(sq_x*he_y + sq_y*he_x) - 70.0f*qd_x*qd_y;
			z[0] = z[0]/c;

			z[1] = y*a - 8.0f*(x*sv_y - y*sv_x) - 56.0f*(qn_x*cu_y - qn_y*cu_x);
			z[1] = z[1]/c;
			break;
		}		
		
		default:
			fprintf(stderr, "unexpected degree\n");
			exit(1);
	}	
}

static void poly_iteration(float x, float y, int *ret, float* roots, int d)
{
	const int max_iters = 128;
	float z[2];
	float zr = x , zi = y;
	float sq_norm, root_sq_norm;
  	int i;
	ret[0] = 0; ret[1] = max_iters;
	for ( i = 0; i < max_iters; i++){
		poly_compute(zr, zi, z, d);
		zr += -z[0];
		zi += -z[1];
		
		if (fabs(zr) > 1e10 || fabs(zi) > 1e10){
			ret[0] = 0;
			ret[1] = i+1;
			break;
		}
		sq_norm = zr*zr + zi*zi;
		
		if (sq_norm < 1e-6){
			ret[0] = 0;
			ret[1] = i+1;
			break;
		}
		
		for (int j = 0; j < d; j++){
			root_sq_norm = SQ(zr - roots[2*j]) + SQ(zi - roots[2*j+1]);
			if (root_sq_norm < 1e-5){
				ret[0] = j+1;
				ret[1] = i+1;
				i = max_iters;
				break;
			}
		}          
	}

}

int thrd_fun(void *args)
{															
  	const thrd_info_t *thrd_info = (thrd_info_t*) args;
  	int **attractors 	= thrd_info->attractors;
  	int **convergences 	= thrd_info->convergences;
  	const int ib 		= thrd_info->ib;
  	const int istep		= thrd_info->istep;
  	const int sz 		= thrd_info->sz;
  	const int d 		= thrd_info->d;
  	const int tx 		= thrd_info->tx;
  	mtx_t *mtx 			= thrd_info->mtx;
  	cnd_t *cnd 			= thrd_info->cnd;
  	int_padded *status 	= thrd_info->status;
	const float dxy   	= 4.0 /sz;
	float x, y;
	int root[2], r_iter, r_index;
	short iters;
	
	float roots[d*2];
	float pi = 3.141596536;
	
	for (int r = 0; r < d; r++)
	{
		roots[2*r] 		= cos(r * 2.0f * pi/d);
		roots[2*r + 1] 	= sin(r * 2.0f * pi/d);	
	}
	
	for ( int ix = ib; ix < sz; ix += istep ) {
		// We allocate the rows of the result before computing, and free them in another thread.
		int *attr_ix 	= (int*) malloc(sz*sizeof(int));
		int *conv_ix 	= (int*) malloc(sz*sizeof(int));
		
		for ( int jx = 0; jx < sz; ++jx ){		
			x = -2.0f + dxy*ix;
			y =  2.0f - dxy*jx;
			poly_iteration(x, y, root, roots, d);
			r_index = root[0];	
			r_iter 	= root[1];
			attr_ix[jx] = r_index;
			conv_ix[jx] = r_iter;
		}

		mtx_lock(mtx);
		attractors[ix] 		= attr_ix; // array of RGB
		convergences[ix] 	= conv_ix; // array of iterations done
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
  int **attractors 		= thrd_info->attractors;
  int **convergences 	= thrd_info->convergences;
  int max_col_val 	  	= thrd_info->max_col_val;
  FILE *file_conv     	= thrd_info->file_conv;
  FILE *file_attr     	= thrd_info->file_attr;
  const int sz 			= thrd_info->sz;
  const int d 			= thrd_info->d;
  const int nthrds 		= thrd_info->nthrds;
  mtx_t *mtx 			= thrd_info->mtx;
  cnd_t *cnd 			= thrd_info->cnd;
  int_padded *status 	= thrd_info->status;

  // create colormaps //
  char clrs_grey[10000];
  clrs_grey[0] = '\0';
  // greyvalues
  for (int ix = 0; ix <= 254; ix += 2) {
      char line[14];
      snprintf(line, sizeof(line), "%i %i %i\n", ix, ix, ix);
      strcat(clrs_grey, line);
  }
  // rgb
  char clrs_rgb[] = "255 100 100\n100 100 255\n100 255 100\n100 255 255\n100 100 100\n255 100 255\n255 255 100\n255 255 255\n100 150 250\n250 150 100\n100 255 200\n";

  write_header(file_conv, sz, max_col_val);
  write_header(file_attr, sz, max_col_val);
	
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

    char *row_str = (char*) malloc(20*sz*sizeof(char));

    // We do not initialize ix in this loop, but in the outer one.
    for ( ; ix < ibnd; ++ix ) {
		write_conv(file_conv, convergences[ix], sz, clrs_grey, row_str);
		write_attr(file_attr, attractors[ix], sz, d, clrs_rgb, row_str);
		free(attractors[ix]);
      	free(convergences[ix]);
    }
    free(row_str);
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
	
	int max_color_val = 255;

	// create filenames
	char f_name_conv[35] = "newton_convergence_x";
	char f_name_attr[35] = "newton_attractors_x";
	strcat(f_name_conv, argv[argc - 1]);
	strcat(f_name_conv, ".ppm");
	strcat(f_name_attr, argv[argc - 1]);
	strcat(f_name_attr, ".ppm");
	//

	// open files
	FILE *file_conv = fopen( f_name_conv, "w" );
	if (file_conv == NULL){
		printf("Error opening file\n");
		return -1;
	}
	FILE *file_attr = fopen( f_name_attr, "w" );
	if (file_attr == NULL){
		printf("Error opening file\n");
		return -1;
	}
	
	assert(n_size != 0 && n_threads != 0);	
	
	int **attractors 	= (int**) malloc(n_size*sizeof(int*));
	int **convergences 	= (int**) malloc(n_size*sizeof(int*));

	thrd_t thrds[n_threads];
	thrd_info_t thrds_info[n_threads];

	thrd_t thrd_check;
	thrd_info_check_t thrd_info_check;
	
	mtx_t mtx;
	mtx_init(&mtx, mtx_plain);

	cnd_t cnd;
	cnd_init(&cnd);

	int_padded status[n_threads];
	int_padded status_wr;

	for ( int tx = 0; tx < n_threads; ++tx ) {
		thrds_info[tx].attractors 	= attractors;
		thrds_info[tx].convergences	= convergences;
		thrds_info[tx].ib	 		= tx;
		thrds_info[tx].istep 		= n_threads;
		thrds_info[tx].sz 			= n_size;
		thrds_info[tx].d 			= n_d;
		thrds_info[tx].tx 			= tx;
		thrds_info[tx].mtx 			= &mtx;
		thrds_info[tx].cnd 			= &cnd;
		thrds_info[tx].status 		= status;
		status[tx].val 				= 0;

		int r = thrd_create(thrds+tx, thrd_fun, (void*) (thrds_info+tx));
		if ( r != thrd_success ) {
			fprintf(stderr, "failed to create thread\n");
			exit(1);
		}
    	thrd_detach(thrds[tx]);
	}	


  	{
    thrd_info_check.attractors 		= attractors;
    thrd_info_check.convergences 	= convergences;
	thrd_info_check.max_col_val 	= 255;
	thrd_info_check.file_conv   	= file_conv;
	thrd_info_check.file_attr   	= file_attr;
    thrd_info_check.sz 				= n_size;
    thrd_info_check.d 				= n_d;
    thrd_info_check.nthrds 			= n_threads;
    thrd_info_check.mtx 			= &mtx;
    thrd_info_check.cnd 			= &cnd;
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

	fclose(file_attr);
	fclose(file_conv);

	free(attractors);
	free(convergences);

	mtx_destroy(&mtx);
	cnd_destroy(&cnd);
	
	return 0;
}
