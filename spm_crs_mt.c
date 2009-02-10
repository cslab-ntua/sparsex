#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>

#include "macros.h"
#include "mt_lib.h"
#include "spm_mt.h"
#include "spm_crs.h"
#include "spm_crs_mt.h"
#include "spmv_method.h"

void *SPM_CRS_MT_NAME(_init_mmf)(char *mmf_file,
                                 unsigned long *rows_nr, unsigned long *cols_nr,
                                 unsigned long *nz_nr)
{
	int i;
	unsigned int nr_cpus, *cpus;
	spm_mt_t *spm_mt;
	spm_mt_thread_t *spm_thread;
	SPM_CRS_MT_TYPE *crs_mt;
	unsigned long cur_row, elems_limit, elems_total=0;

	SPM_CRS_TYPE *crs;
	crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);
	mt_get_options(&nr_cpus, &cpus);

	spm_mt = malloc(sizeof(spm_mt_t));
	if ( !spm_mt ){
		perror("malloc");
		exit(1);
	}

	spm_mt->nr_threads = nr_cpus;
	spm_mt->spm_threads = malloc(sizeof(spm_mt_thread_t)*nr_cpus);
	if ( !spm_mt->spm_threads ){
		perror("malloc");
		exit(1);
	}

	crs_mt = malloc(sizeof(SPM_CRS_MT_TYPE)*nr_cpus);
	if ( !crs_mt ){
		perror("malloc");
		exit(1);
	}

	for (i=0, cur_row=0; i<nr_cpus; i++){
		elems_limit = (*nz_nr - elems_total) / (nr_cpus - i);
		spm_thread = spm_mt->spm_threads+ i;
		spm_thread->cpu = cpus[i];
		spm_thread->spm = crs_mt + i;

		unsigned long elems, rows;
		crs_mt[i].row_start = cur_row;
		for (elems=0,rows=0; ; ){
			elems += crs->row_ptr[cur_row+1] - crs->row_ptr[cur_row];
			cur_row++;
			//printf("i:%d nr_cpus:%d cur_row:%lu rows_nr:%lu elems:%lu elems_limit:%lu\n", i, nr_cpus, cur_row, *rows_nr, elems, elems_limit);
			#if 0
			if (i != (nr_cpus -1)){
				if ( elems >= elems_limit )
					break;
			} else {
				if (cur_row == *rows_nr)
					break;
			}
			#endif
			if (elems >= elems_limit){
				break;
			}
		}

		elems_total += elems;
		crs_mt[i].row_end = cur_row;
		crs_mt[i].nnz_nr = elems;
		crs_mt[i].crs = crs;
	}

	free(cpus);
	return spm_mt;
}

void SPM_CRS_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	SPM_CRS_IDX_TYPE *row_ptr = crs_mt->crs->row_ptr;
	SPM_CRS_IDX_TYPE *col_ind = crs_mt->crs->col_ind;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	register ELEM_TYPE yr;

	unsigned long i,j;
	for (i=row_start; i<row_end; i++){
		yr = (ELEM_TYPE)0;
		for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
			yr += (values[j] * x[col_ind[j]]);
		}
		y[i] = yr;
	}
}
XSPMV_MT_METH_INIT(
 SPM_CRS_MT_NAME(_multiply),
 SPM_CRS_MT_NAME(_init_mmf),
 NULL,
 NULL,
 sizeof(ELEM_TYPE)
)
