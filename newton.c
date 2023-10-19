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

#define MIN(a,b) ({ __typeof__ (a) _a = (a); \
										__typeof__ (b) _b = (b); \
										_a < _b ? _a : _b; })

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

typedef struct {
  int **w;
	int **f;
  int max_col_val;
  FILE *file_conv;
  FILE *file_attr;
  int ib;
  int istep;
  int sz;
	int d;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_wr_info_t;

static void
write_header(FILE *file, int n_size, int max_color_val){
    fprintf(file, "%s\n%i %i\n%i\n", "P3", n_size, n_size, max_color_val);
}

static void
write_conv(FILE *file, int* convergence, int n_size)
{
  int color_str_len = 20;
  int row_str_len_sum = 0;
  int offset = 0;
  char *row_str = (char*) malloc(n_size*color_str_len*sizeof(char));
  char colors[10000];
  colors[0] = '\0';

  // create colormap for greyvalues
  for (int ix = 0; ix <= 254; ix += 2) {
      char line[14];
      snprintf(line, sizeof(line), "%i %i %i\n", ix, ix, ix);
      strcat(colors, line);
  }

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

  fwrite( row_str, sizeof(char), row_str_len_sum, file);
  free(row_str);
  free(convergence);
}

static void
write_attr(FILE *file, int* attractor, int n_size, int n_degree)
{
    int color_str_len = 12;

    // color map
    char c_0[] = "100 100 100\n";
    char c_1[] = "100 100 255\n";
    char c_2[] = "100 255 100\n";
    char c_3[] = "100 255 255\n";
    char c_4[] = "255 100 100\n";
    char c_5[] = "255 100 255\n";
    char c_6[] = "255 255 100\n";
    char c_7[] = "255 255 255\n";
    char c_8[] = "100 150 250\n";
    char c_9[] = "250 150 100\n";

    char colors[140];
    sprintf(colors, "%s%s%s%s%s%s%s%s%s%s", c_0, c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9);

    char *row_str = (char*) malloc(n_size*color_str_len*sizeof(char));

    for ( size_t ix = 0, jx = 0; jx < n_size; ix += color_str_len, ++jx ){
      assert( attractor[jx] <= 9 && attractor[jx] >= 0);
      memcpy( row_str + ix, colors + attractor[jx]*color_str_len, color_str_len);        
    }
    
    fwrite( row_str, sizeof(char), n_size*color_str_len, file);
    free(row_str);
}

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

static void poly_iteration(float x, float y, int *ret, int d)
{
	float roots_x[d];
	float roots_y[d];
	float pi = 3.141596536;
	const int max_iters = 100;
	for (int r = 0; r < d; r++)
	{
		roots_x[r] = cos(r * 2.0f * pi/d);
		roots_y[r] = sin(r * 2.0f * pi/d);	
	}
	float z[2];
	float zr = x , zi = y;
	float sq_norm, root_sq_norm;
	for (int i = 0; i < max_iters; i++){
		poly_compute(zr, zi, z, d);
		zr += -z[0];
		zi += -z[1];
		
		if (fabs(zr) > 1e10 || fabs(zi) > 1e10){
			ret[0] = 0;
			ret[1] = i;
			break;
		}
		sq_norm = zr*zr + zi*zi;
		
		if (sq_norm < 1e-6){
			ret[0] = 0;
			ret[1] = i;
			break;
		}
		
		for (int j = 0; j < d; j++){
			root_sq_norm = SQ(zr - roots_x[j]) + SQ(zi - roots_y[j]);
			if (root_sq_norm < 1e-6){
				ret[0] = j+1;
				ret[1] = i;
				i = max_iters;
				break;
			}
		}
	}	
}

int thrd_write(void *args)
{
  const thrd_wr_info_t *thrd_wr_info = (thrd_wr_info_t*) args;
  int **w 					  = thrd_wr_info->w;
  int **f 						= thrd_wr_info->f;
  int max_col_val 	  = thrd_wr_info->max_col_val;
  FILE *file_conv     = thrd_wr_info->file_conv;
  FILE *file_attr     = thrd_wr_info->file_attr;
  const int ib 				= thrd_wr_info->ib;
  const int istep			= thrd_wr_info->istep;
  const int sz 				= thrd_wr_info->sz;
  const int d 				= thrd_wr_info->d;
  const int tx 				= thrd_wr_info->tx;
  mtx_t *mtx 					= thrd_wr_info->mtx;
  cnd_t *cnd 					= thrd_wr_info->cnd;
  int_padded *status 	= thrd_wr_info->status;

  write_header(file_conv, sz, max_col_val);
  write_header(file_attr, sz, max_col_val);

  for ( int ix = 0; ix < sz; ++ix ){
    write_conv(file_conv, f[ix], sz);
  }

  for ( int ix = 0; ix < sz; ++ix ){
    write_attr(file_attr, w[ix], sz, d);
  }
  
  return 0;
}

int thrd_fun(void *args)
{
	const float rgb_colors[33] = {255., 255., 255.,
																0., 0., 255., 		255., 255., 0,\
																255., 0., 255., 	0., 128., 128.,\
																51.f, 153., 102, 	204., 250., 180.,\
																153., 51., 0., 		120., 30., 200.,\
																255., 0., 0.,			255., 255., 200.};
																
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
	float x, y;
	int root[2], r_iter, r_index;
	short iters;
	for ( int ix = ib; ix < sz; ix += istep ) {
		// We allocate the rows of the result before computing, and free them in another thread.
		float *wix 	= (float*) malloc(3*sz*sizeof(float));
		int *fix 		= (int*) 	malloc(sz*sizeof(int));
		
		for ( int jx = 0; jx < sz; ++jx ){		
			x = -2.0f + dxy*ix;
			y =  2.0f - dxy*jx;
			poly_iteration(x, y, root, d);
			r_index = root[0];	
			r_iter 	= root[1];
			wix[3 * jx] 			= rgb_colors[3 * r_index]; 	// R
			wix[3 * jx + 1] 	= rgb_colors[3 * r_index + 1]; // G
			wix[3 * jx + 2] 	= rgb_colors[3 * r_index + 2]; // B
			fix[jx] 					= MIN(iters * 255.0f/100.0f, 255.0f);
			printf("%d\n", r_iter);
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
      // free(f[ix]);
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

  thrd_t thrd_wr;
  thrd_wr_info_t thrd_wr_info;

  thrd_t thrd_check;
  thrd_info_check_t thrd_info_check;
  
  mtx_t mtx;
  mtx_init(&mtx, mtx_plain);

  cnd_t cnd;
  cnd_init(&cnd);

  int_padded status[n_threads];
  int_padded status_wr;

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

	int max_color_val = 255;

  // create filename
  char f_name_conv[] = "newton_convergence_x";
  char f_name_attr[] = "newton_attractors_x";
  strcat(f_name_conv, argv[argc - 1]);
  strcat(f_name_conv, ".ppm");
  strcat(f_name_attr, argv[argc - 1]);
  strcat(f_name_attr, ".ppm");
  //

  // open file
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

  // create temporary attractors
  int **attractors  = (int**) malloc(n_size*sizeof(int*));
  int *attractor    = (int*) malloc(n_size*n_size*sizeof(int));

  for ( int ix = 0; ix < n_size*n_size; ++ix)
    attractor[ix] = (int) rand() % 10;

  for ( int ix = 0; ix < n_size; ++ix)
    attractors[ix] = attractor + n_size*ix;
  
  {
  thrd_wr_info.w 			      = attractors;
  thrd_wr_info.f 			      = f;
  thrd_wr_info.max_col_val  = 255;
  thrd_wr_info.file_conv    = file_conv;
  thrd_wr_info.file_attr    = file_attr;
  thrd_wr_info.ib	 		      = 0;
  thrd_wr_info.istep 	      = n_threads;
  thrd_wr_info.sz 		      = n_size;
  thrd_wr_info.d 			      = n_d;
  thrd_wr_info.tx 		      = 0;
  thrd_wr_info.mtx 		      = &mtx;
  thrd_wr_info.cnd 		      = &cnd;
  thrd_wr_info.status       = status;
  status_wr.val 	          = 0;

  int r = thrd_create(&thrd_wr, thrd_write, (void*) (&thrd_wr_info));
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }

  {
    int r;
    thrd_join(thrd_wr, &r);
  }

  
  free(attractor);
  free(attractors);

  fclose(file_attr);
  fclose(file_conv);

  free(w);
  free(f);

  mtx_destroy(&mtx);
  cnd_destroy(&cnd);
	
	return 0;
}
