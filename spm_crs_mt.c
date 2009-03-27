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

void SPM_CRS_MT_NAME(_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
	SPM_CRS_NAME(_destroy)(crs_mt->crs);
	free(crs_mt);
	free(spm_thread);
	free(spm_mt);
}

uint64_t SPM_CRS_MT_NAME(_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
	return SPM_CRS_NAME(_size)(crs_mt->crs);
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
 SPM_CRS_MT_NAME(_size),
 SPM_CRS_MT_NAME(_destroy),
 sizeof(ELEM_TYPE)
)

#ifdef SPM_NUMA

#include <numa.h>

static int numa_node_from_cpu(int cpu)
{
	struct bitmask *bmp;
	int nnodes, node, ret;

	bmp = numa_allocate_cpumask();
	nnodes =  numa_num_configured_nodes();
	for (node = 0; node < nnodes; node++){
		numa_node_to_cpus(node, bmp);
		if (numa_bitmask_isbitset(bmp, cpu)){
			ret = node;
			goto end;
		}
	}
	ret = -1;
end:
	numa_bitmask_free(bmp);
	return ret;
}

void *SPM_CRS_MT_NAME(_numa_init_mmf)(char *mmf_file,
                                      uint64_t *rows_nr, uint64_t *cols_nr,
                                      uint64_t *nz_nr)
{
	spm_mt_t *spm_mt;
	int i;
	spm_mt = SPM_CRS_MT_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);

	/* sanity check */
	if (numa_available() == -1){
		perror("numa not available");
		exit(1);
	}

	/* keep a reference to original crs */
	SPM_CRS_TYPE *crs = ((SPM_CRS_MT_TYPE *)spm_mt->spm_threads->spm)->crs;
	ELEM_TYPE *values = crs->values;
	SPM_CRS_IDX_TYPE *colind = crs->col_ind;
	SPM_CRS_IDX_TYPE *rowptr = crs->row_ptr;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
		/* get numa node from cpu */
		int node = numa_node_from_cpu(spm_thread->cpu);
		/* allocate new space */
		SPM_CRS_TYPE *numa_crs = numa_alloc_onnode(sizeof(SPM_CRS_TYPE), node);
		if (!numa_crs){
			perror("numa_alloc_onnode");
			exit(1);
		}
		uint64_t nnz = crs_mt->nnz_nr;
		uint64_t nrows = crs_mt->row_end - crs_mt->row_start;
		numa_crs->values = numa_alloc_onnode(sizeof(ELEM_TYPE)*nnz, node);
		numa_crs->col_ind = numa_alloc_onnode(sizeof(SPM_CRS_IDX_TYPE)*nnz, node);
		numa_crs->row_ptr = numa_alloc_onnode(sizeof(SPM_CRS_IDX_TYPE)*(nrows+1), node);
		if (!numa_crs->values || !numa_crs->col_ind || !numa_crs->row_ptr){
			perror("numa_alloc_onnode");
			exit(1);
		}
		/* copy data */
		memcpy(numa_crs->values, values, sizeof(ELEM_TYPE)*nnz);
		values += nnz;
		memcpy(numa_crs->col_ind, colind, sizeof(SPM_CRS_IDX_TYPE)*nnz);
		colind += nnz;
		memcpy(numa_crs->row_ptr, rowptr, sizeof(SPM_CRS_IDX_TYPE)*(nrows+1));
		rowptr += nrows;
		/* make the swap */
		numa_crs->nz = nnz;
		numa_crs->nrows = nrows;
		numa_crs->ncols = crs->ncols;
		crs_mt->crs = numa_crs;
	}
	SPM_CRS_NAME(_destroy)(crs);

	return spm_mt;
}

void SPM_CRS_MT_NAME(_numa_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	int i;
	SPM_CRS_MT_TYPE *crs_mt0 = (SPM_CRS_MT_TYPE *)spm_mt->spm_threads->spm;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
		SPM_CRS_TYPE *numa_crs = crs_mt->crs;
		uint64_t nnz = numa_crs->nz;
		uint64_t nrows = crs_mt->row_end - crs_mt->row_start;
		numa_free(numa_crs->values, sizeof(ELEM_TYPE)*nnz);
		numa_free(numa_crs->col_ind, sizeof(SPM_CRS_IDX_TYPE)*nnz);
		numa_free(numa_crs->row_ptr, sizeof(SPM_CRS_IDX_TYPE)*(nrows+1));
		numa_free(numa_crs, sizeof(SPM_CRS_TYPE));
	}
	free(spm_mt->spm_threads);
	free(crs_mt0);
	free(spm_mt);
}

uint64_t SPM_CRS_MT_NAME(_numa_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	int i;
	uint64_t ret = 0;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread  = spm_mt->spm_threads + i;
		SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
		ret += SPM_CRS_NAME(_size)(crs_mt->crs);

	}
	return ret;
}

void SPM_CRS_MT_NAME(_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	SPM_CRS_IDX_TYPE *row_ptr = crs_mt->crs->row_ptr;
	SPM_CRS_IDX_TYPE *col_ind = crs_mt->crs->col_ind;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	const unsigned long nrows = row_end - row_start;
	register ELEM_TYPE yr;
	ELEM_TYPE *y = out->elements + row_start;

	unsigned long i,j;
	uint64_t j0 = row_ptr[0];
	for (i=0; i<nrows; i++){
		yr = (ELEM_TYPE)0;
		for (j=row_ptr[i]-j0; j<row_ptr[i+1]-j0; j++) {
			yr += ( *(values++) * x[col_ind[j]]);
		}
		y[i] = yr;
	}
}

XSPMV_MT_METH_INIT(
 SPM_CRS_MT_NAME(_numa_multiply),
 SPM_CRS_MT_NAME(_numa_init_mmf),
 SPM_CRS_MT_NAME(_numa_size),
 SPM_CRS_MT_NAME(_numa_destroy),
 sizeof(ELEM_TYPE)
)

#endif /* SPM_NUMA */

