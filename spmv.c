#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <math.h>

#include "tsc.h"
#include "prfcnt.h"
#include "matrix.h"
#include "method.h"
#include "spm_crs.h"
//#include "prefetch.h"

#ifndef LOOPS
#define LOOPS 128
#endif
static unsigned long ncols;

static void do_setaffinity(cpu_set_t *mask)
{
	int err;

	err = sched_setaffinity(0, sizeof(cpu_set_t), mask);
	if (err){
		perror("sched_setaffinity");
		exit(1);
	}
}


static inline int elems_neq(elem_t a, elem_t b)
{
	if ( fabs((double)(a-b)/(double)a)  > 0.0000001 ){
		return 1;
	}
	return 0;
}

typedef void spmv_fn_t(void *matrix, vector_t *in, vector_t *out);
typedef void *spmv_load_fn_t(char *mmf_file, 
                             unsigned long *nrows, unsigned long *ncols, 
			     unsigned long *nnz);

int check_spm_crs(char *mmf_file, vector_t *x, vector_t *y)
{
	vector_t *y_correct;
	unsigned long i;
	method_t  *m;
	spmv_fn_t *spmv_fn;
	spmv_load_fn_t *load_fn;
	unsigned long rows_nr, cols_nr, nz_nr;
	spm_crs64_t *matrix;

	printf("checking ...\n");
	m = method_get("spm_crs64_multiply");
	spmv_fn = m->fn;
	load_fn = m->data;
	matrix = load_fn(mmf_file, &rows_nr, &cols_nr, &nz_nr);
	y_correct = vector_create(cols_nr);
	vector_init(y_correct, (elem_t)0);
	spmv_fn(matrix, x, y_correct);

	if ( y_correct->size != y->size ){
		fprintf(stderr, "check_spm_crs: sizes differ\n");
		exit(1);
	}

	for (i=0; i < y_correct->size; i++){
		if ( elems_neq(y_correct->elements[i], y->elements[i]) ){
			fprintf(stderr, "element %ld differs: %10.20lf != %10.20lf\n", 
			                 i,  y_correct->elements[i], y->elements[i]);
			exit(1);
		}
	}

	return 0;
}

int main(int argc, char **argv)
{
	int fd, i, ret;
	int cpu=0;
	vector_t *x,*y;
	spmv_fn_t *fn;
	spmv_load_fn_t *load_fn;
	method_t *m;
	char *method;
	char *ncols_str;
	unsigned long rows_nr, cols_nr, nz_nr;
	tsc_t tsc;
	prfcnt_t prfcnt;
	cpu_set_t cpu_mask;
	void *matrix;

	ncols_str = getenv("SPMV_NCOLS"); 
	if ( ncols_str ){
		ncols = atol(ncols_str);
	} else {
		ncols = ~0; // kill dah bitch
	}

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmarket_file [method]\n", argv[0]);
		method_fprint(stderr, "available methods:\n", "\t", "\n", "");
		exit(1);
	}

	method = (argc > 2) ? argv[2] : "spm_crs_multiply_safe";
	m = method_get(method);
	if ( !m ) {
		fprintf(stderr, "No such function: %s\n", argv[2]);
		method_fprint(stderr, "available methods:\n", "\t", "\n", "");
		exit(1);
	}
	fn = m->fn;
	load_fn = m->data;

	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu, &cpu_mask);
	do_setaffinity(&cpu_mask);

	matrix = load_fn(argv[1], &rows_nr, &cols_nr, &nz_nr);
	x = vector_create(cols_nr);
	y = vector_create(cols_nr);
	
	srand(time(NULL));
	tsc_init(&tsc);
	prfcnt_init(&prfcnt,cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);
	for (i=0; i<LOOPS ;  i++){
		vector_init_rand_range(x, (elem_t)-1000, (elem_t)1000);
		//vector_init(x, (elem_t)3);
		vector_init(y, (elem_t)0);
		//prfcnt_start(&prfcnt);
		tsc_start(&tsc);
		fn(matrix, x, y);
		tsc_pause(&tsc);
		//prfcnt_pause(&prfcnt);
		#ifdef SPMV_CHECK
		//vector_print(x);
		check_spm_crs(argv[1], x, y);
		#endif
	}
	//tsc_report(&tsc);
	tsc_shut(&tsc);
	const float secs = tsc_getsecs(&tsc);
	const float mf = (((float)(nz_nr*2*LOOPS))/secs)/(1000*1000);
	const float mb = (((float)(nz_nr*2+rows_nr*2+ cols_nr)*8*LOOPS)/secs)/(1024.0*1024.0);
	printf("%s %s %f %f %f\n", basename(argv[1]), method, mf, mb, secs);
	//prfcnt_report(&prfcnt);

	return 0;
}

void do_check_loop(char *file, char *method, void *matrix,
                   unsigned long loops, unsigned long cols_nr)
{
	unsigned long i;
	method_t *m;
	vector_t *x,*y;
	spmv_fn_t *fn;

	x = vector_create(cols_nr);
	y = vector_create(cols_nr);

	m = method_get(method);
	if ( !m ) {
		fprintf(stderr, "No such function: %s\n", method);
		method_fprint(stderr, "available methods:\n", "\t", "\n", "");
		exit(1);
	}
	fn = m->fn;
	for (i=0; i<loops ;  i++){
		vector_init_rand_range(x, (elem_t)-1000, (elem_t)1000);
		//vector_init(x, (elem_t)3);
		vector_init(y, (elem_t)0);
		fn(matrix, x, y);
		//vector_print(x);
		check_spm_crs(file, x, y);
	}

	vector_destroy(x);
	vector_destroy(y);
}

void do_bench_loop(char *file, char *method, void *matrix,
		  unsigned long loops,
                  unsigned long cols_nr, unsigned long nz_nr)
{
	unsigned long i;
	method_t *m;
	vector_t *x,*y;
	spmv_fn_t *fn;
	tsc_t tsc;

	x = vector_create(cols_nr);
	y = vector_create(cols_nr);

	m = method_get(method);
	if ( !m ) {
		fprintf(stderr, "No such function: %s\n", method);
		method_fprint(stderr, "available methods:\n", "\t", "\n", "");
		exit(1);
	}
	fn = m->fn;

	tsc_init(&tsc);
	for (i=0; i<loops; i++){
		vector_init_rand_range(x, (elem_t)-1000, (elem_t)1000);
		vector_init(y, (elem_t)0);
		tsc_start(&tsc);
		fn(matrix, x, y);
		tsc_pause(&tsc);
	}
	tsc_shut(&tsc);
	const float secs = tsc_getsecs(&tsc);
	const float mf = (((float)(nz_nr*2*LOOPS))/secs)/(1000*1000);
	//const float mb = (((float)(nz_nr*2+rows_nr*2+cols_nr)*8*LOOPS)/secs)/(1024.0*1024.0);
	printf("%s %s %f %f\n", basename(file), method, mf, secs);
}

#if 0
void spm_crs_multiply_safe(spm_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->crs.values;
	index_t *row_ptr = matrix->crs.row_ptr;
	index_t *col_ind = matrix->crs.col_ind;
	unsigned long n = matrix->crs.nrows;

	unsigned long i, j;
	for(i=0; i<n; i++) {
		//printf("row: %lu (%lu,%lu)\n", i, row_ptr[i], row_ptr[i+1]);
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) { 
			y[i] += (values[j] * x[col_ind[j]]);
		}
	}
}
METHOD_INIT(spm_crs_multiply_safe, spm_crs_load)

void spm_crs_multiply_nocomp(spm_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->crs.values;
	index_t *row_ptr = matrix->crs.row_ptr;
	index_t *col_ind = matrix->crs.col_ind;
	unsigned long n = matrix->crs.nrows;

	unsigned long i, j;
	for(i=0; i<n; i++) {
		//__preloadq(&y[i]);
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) { 
			//__preloadq(&values[j]);
			//__preloadq(&x[col_ind[j]]);
		}
	}
}
METHOD_INIT(spm_crs_multiply_nocomp, spm_crs_load)

void spm_crs_multiply_nocomp_v(spm_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->crs.values;
	index_t *row_ptr = matrix->crs.row_ptr;
	index_t *col_ind = matrix->crs.col_ind;
	unsigned long n = matrix->crs.nrows;

	unsigned long i, j;
	for(i=0; i<n; i++) {
		//__preloadq(&y[i]);
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) { 
			//__preloadq(&values[j]);
		}
	}
}
METHOD_INIT(spm_crs_multiply_nocomp_v, spm_crs_load)

void spm_crs_multiply_nocomp_x(spm_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->crs.values;
	index_t *row_ptr = matrix->crs.row_ptr;
	index_t *col_ind = matrix->crs.col_ind;
	unsigned long n = matrix->crs.nrows;

	unsigned long i, j;
	for(i=0; i<n; i++) {
		//__preloadq(&y[i]);
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) { 
			//__preloadq(&x[col_ind[j]]);
		}
	}
}
METHOD_INIT(spm_crs_multiply_nocomp_x, spm_crs_load)

static void spm_crs_multiply_no_indrct1(spm_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->crs.values;
	index_t *row_ptr = matrix->crs.row_ptr;
	index_t *col_ind = matrix->crs.col_ind;
	unsigned long n = matrix->crs.nrows;

	unsigned long i, j, j_max = 0, j_prev=0;
	printf("ncols=%lu\n", ncols);

	for(i=0; i<n; i++) {
		for(j=ncols*i; j < ncols*(i+1); j++) { 
			y[i] += (values[j] * x[col_ind[j]]);
		}
	}
}
METHOD_INIT(spm_crs_multiply_no_indrct1, spm_crs_load)

static void spm_crs_multiply_no_indrct2(spm_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->crs.values;
	index_t *row_ptr = matrix->crs.row_ptr;
	index_t *col_ind = matrix->crs.col_ind;
	unsigned long n = matrix->crs.nrows;
	unsigned long ncols = matrix->crs.ncols;

	unsigned long i, j, k, j_max = 0, j_prev=0;
	//ncols = 1;

	for(i=0; i<n; i++) {
		for(j=row_ptr[i], k=0; j < row_ptr[i+1]; j++, k++) { 
			y[i] += (values[j] * x[k]);
		}
	}
}
METHOD_INIT(spm_crs_multiply_no_indrct2, spm_crs_load)
#endif
