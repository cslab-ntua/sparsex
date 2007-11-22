#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>

#include "tsc.h"
#include "prfcnt.h"
#include "matrix.h"
#include "sparse.h"

#define LOOPS 128

static void do_setaffinity(cpu_set_t *mask)
{
	int err;

	err = sched_setaffinity(0, sizeof(cpu_set_t), mask);
	if (err){
		perror("sched_setaffinity");
		exit(1);
	}
}

typedef void spmv_fn_t(spm_crs_t *matrix, vector_t *in, vector_t *out);

static void spm_crs_multiply_1(spm_crs_t *matrix, vector_t *in, vector_t *out);

static struct {
	char      *name;
	spmv_fn_t *fn;
} spmv_methods[] = {
	{"crs_spmv", spm_crs_multiply},
	{"crs_spmv1", spm_crs_multiply_1},
	{NULL, NULL}
};


static spmv_fn_t *get_method(char *name)
{
	int i;
	spmv_fn_t *ret = NULL;
	for ( i=0; spmv_methods[i].name; i++){
		if ( strcmp(spmv_methods[i].name, name) == 0){
			ret = spmv_methods[i].fn;
			break;
		} 
	}

	return ret;
}

int main(int argc, char **argv)
{
	int fd, i, ret;
	int cpu=0;
	spm_crs_t *crs;
	vector_t *x,*y;
	spmv_fn_t *fn;
	char *method;
	unsigned long rows_nr, cols_nr, nz_nr;
	tsc_t tsc;
	prfcnt_t prfcnt;
	cpu_set_t cpu_mask;

	if (argc < 2){
		fprintf(stderr, "Usage: %s mmarket_file [method]\n", argv[0]);
		exit(1);
	}

	method = (argc > 2) ? argv[2] : "crs_spmv";
	fn = get_method(method);
	if ( !fn ) {
		fprintf(stderr, "No such function: %s\n", argv[2]);
		exit(1);
	}

	fd = spm_open_mmarket_file(argv[1]);
	ret = spm_read_header(fd, &rows_nr, &cols_nr, &nz_nr);
	if ( ret < 0){
		fprintf(stderr, "Can't read header\n");
		exit(1);
	}

	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu, &cpu_mask);
	do_setaffinity(&cpu_mask);

	crs = spm_crs_create(rows_nr, cols_nr, nz_nr);
	spm_crs_init(crs, fd);
	spm_crs_mkdummy(crs, SPM_DUMMY_COLIND);
	x = vector_create(cols_nr);
	y = vector_create(cols_nr);
	
	srand(time(NULL));
	tsc_init(&tsc);
	prfcnt_init(&prfcnt,cpu,PRFCNT_FL_T0|PRFCNT_FL_T1);
	for (i=0; i<LOOPS ;  i++){
		vector_init_rand_range(x, (elem_t)-1000, (elem_t)1000);
		vector_init(y, (elem_t)0);
		prfcnt_start(&prfcnt);
		tsc_start(&tsc);
		fn(crs, x, y);
		tsc_pause(&tsc);
		prfcnt_pause(&prfcnt);
	}
	//tsc_report(&tsc);
	tsc_shut(&tsc);
	float mflops = (((float)(nz_nr*2*LOOPS))/tsc_getsecs(&tsc))/(1000*1000);
	printf("%s %s %f\n", basename(argv[1]), method, mflops);
	prfcnt_report(&prfcnt);

	return 0;
}

static void spm_crs_multiply_1(spm_crs_t *matrix, vector_t *in, vector_t *out)
{
	elem_t *y = out->elements;
	elem_t *x = in->elements;
	elem_t *values = matrix->values;
	index_t *row_ptr = matrix->row_ptr;
	index_t *col_ind = matrix->col_ind;
	unsigned long n = matrix->nrows;

	unsigned long i, j;

	for(i=0; i<n; i++) {
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) { 
			y[i] += (values[j] * x[col_ind[j]]);
		}
	}
}

