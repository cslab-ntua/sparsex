#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>

#include "tsc.h"
//#include "prfcnt.h"
#include "matrix.h"
#include "dmv_vec.h"

#define LOOPS 128
#define MAT_SIZE 1024


int pollute_cache()
{
	char *b;
	int  i, sum;
	const int size = (CACHE_BYTES << 1);

	b = malloc(size);
	if (!b){
		perror("malloc");
		exit(1);
	}

	for (i=0 ; i < size; i++)
		b[i] = i;

	for (sum=0, i=0 ; i < size; i++)
		sum += b[i];

	free(b);

	return sum;
}

static void do_setaffinity(cpu_set_t *mask)
{
	int err;

	err = sched_setaffinity(0, sizeof(cpu_set_t), mask);
	if (err){
		perror("sched_setaffinity");
		exit(1);
	}
}

static void dmxv_mul(matrix_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	unsigned long nrows = matrix->nrows;
	unsigned long ncols = matrix->ncols;

	unsigned long i, j;

	for (i=0; i<nrows; i++){
		elem_t *v = matrix->rows[i];
		for (j=0; j < ncols; j++){
			y[i] += v[j]*x[j];
		}
	}

}

static void dmxv_mul_vec(matrix_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	unsigned long nrows = matrix->nrows;
	unsigned long ncols = matrix->ncols;

	unsigned long i, j;

	for (i=0; i<nrows; i++){
		elem_t *v = matrix->rows[i];
		y[i] = dmv_vec_va(v, x, ncols);
	}

}

#define MIN(x,y) (x < y) ? x : y
static void dmxv_mul_bvec(matrix_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	unsigned long nrows = matrix->nrows;
	unsigned long ncols = matrix->ncols;

	unsigned long i, j, col_max, ci_start, size;

	size = CACHE_BYTES / (sizeof(elem_t)*2);
	for (ci_start=0; ci_start<ncols; ci_start+=size){
		col_max = MIN(ci_start+size, ncols);
		for ( i=0; i<nrows; i++){
			elem_t *v = &matrix->rows[i][ci_start];
			elem_t *x = &in->elements[ci_start];
			y[i] = dmv_vec(v, x, col_max);
			
		}
	}
}

int main(int argc, char **argv)
{
	int fd, i, ret;
	int cpu=0;
	vector_t *x,*y;
	matrix_t *m;
	tsc_t tsc;
	//prfcnt_t prfcnt;
	cpu_set_t cpu_mask;
	unsigned long nr_rows, nr_cols;

	nr_rows = (argc<2) ? MAT_SIZE: atoi(argv[1]);
	nr_cols = (argc<3) ? MAT_SIZE: atoi(argv[2]);

	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu, &cpu_mask);
	do_setaffinity(&cpu_mask);

	x = vector_create(nr_rows);
	y = vector_create(nr_rows);
	m = matrix_create(nr_rows,nr_cols);
	matrix_init_rand_range(m, (elem_t)-1000, (elem_t)1000);
	
	srand(time(NULL));
	tsc_init(&tsc);
	//prfcnt_init(&prfcnt,cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);
	for (i=0; i<LOOPS ;  i++){
		vector_init_rand_range(x, (elem_t)-1000, (elem_t)1000);
		vector_init(y, (elem_t)0);
		//pollute_cache();
		//prfcnt_start(&prfcnt);
		tsc_start(&tsc);
		dmxv_mul(m, x, y);
		tsc_pause(&tsc);
		//prfcnt_pause(&prfcnt);
	}
	//tsc_report(&tsc);
	tsc_shut(&tsc);
	float mflops = (((float)(nr_rows*nr_cols*2*LOOPS))/tsc_getsecs(&tsc))/(1000*1000);
	printf("%lux%lu %f\n", nr_rows, nr_cols, mflops);
	//prfcnt_report(&prfcnt);

	return 0;
}


