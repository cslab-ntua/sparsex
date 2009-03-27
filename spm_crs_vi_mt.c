#include <stdlib.h>
#include <inttypes.h>

#include "spm_mt.h"
#include "spm_crs_vi.h"
#include "spm_crs_vi_mt.h"
#include "mt_lib.h"

#include "pthread.h"

#define SPM_CRSVI_MT_NAME(name) \
        CON8(spm_crs, SPM_CRSVI_CI_BITS, _vi, SPM_CRSVI_VI_BITS, _, ELEM_TYPE, _mt, name)
#define SPM_CRSVI_MT_TYPE SPM_CRSVI_MT_NAME(_t)

void *SPM_CRSVI_MT_NAME(_init_mmf)(char *mmf_file,
                                   unsigned long *rows_nr, unsigned long *cols_nr,
                                   unsigned long *nz_nr)
{
	unsigned int nr_cpus, *cpus_affinity;
	SPM_CRSVI_TYPE *crsvi;
	SPM_CRSVI_MT_TYPE *crsvi_mt;
	int i;
	unsigned long cur_row, elems_limit, elems_total=0;

	crsvi = SPM_CRSVI_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);
	mt_get_options(&nr_cpus, &cpus_affinity);

	spm_mt_t *spm_mt;
	spm_mt_thread_t *spm_thread;

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

	crsvi_mt = malloc(sizeof(SPM_CRSVI_MT_TYPE)*nr_cpus);
	if ( !crsvi_mt ){
		perror("malloc");
		exit(1);
	}


	for (i=0, cur_row=0; i<nr_cpus; i++){
		elems_limit = (*nz_nr - elems_total) / (nr_cpus - i);
		spm_thread = spm_mt->spm_threads + i;
		spm_thread->cpu = cpus_affinity[i];
		spm_thread->spm = crsvi_mt + i;

		unsigned long elems, rows;
		crsvi_mt[i].row_start = cur_row;
		for (elems=0,rows=0; ; ){
			elems += crsvi->row_ptr[cur_row+1] - crsvi->row_ptr[cur_row];
			cur_row++;
			//printf("i:%d nr_cpus:%d cur_row:%lu rows_nr:%lu elems:%lu elems_limit:%lu\n", i, nr_cpus, cur_row, *rows_nr, elems, elems_limit);
			if (elems >= elems_limit){
				break;
			}
		}

		elems_total += elems;
		crsvi_mt[i].row_end = cur_row;
		crsvi_mt[i].nnz = elems;
		crsvi_mt[i].crsvi = crsvi;
	}

	free(cpus_affinity);

	return spm_mt;
}

void SPM_CRSVI_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crsvi_mt->crsvi->values;
	const register SPM_CRSVI_CI_TYPE *row_ptr = crsvi_mt->crsvi->row_ptr;
	const register SPM_CRSVI_CI_TYPE *col_ind = crsvi_mt->crsvi->col_ind;
	const register SPM_CRSVI_VI_TYPE *val_ind = crsvi_mt->crsvi->val_ind;
	const unsigned long row_start = crsvi_mt->row_start;
	const unsigned long row_end = crsvi_mt->row_end;
	register ELEM_TYPE yr;

	//printf("(%lu) crsvi_mt:%p row_start: %lu row_end: %lu\n", pthread_self(), crsvi_mt, row_start, row_end);
	register unsigned long i, j=0;
	for (i=row_start; i<row_end; i++){
		yr = (ELEM_TYPE)0;
		__asm__ __volatile__ ("# loop start\n\t");
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) {
			yr += (values[val_ind[j]] * x[col_ind[j]]);
			//printf("(%lu) i:%lu j:%lu val_ind:%ld col_ind:%lu v:%lf x:%lf yr:%lf\n", (unsigned long)pthread_self(),i, j, (long)val_ind[j], (unsigned long)col_ind[j], values[val_ind[j]], x[col_ind[j]], yr);
		}
		y[i] = yr;
		//printf("(%lu): y[%lu] = %lf\n", (unsigned long)pthread_self(), i, (double)yr);
		__asm__ __volatile__ ("# loop end\n\t");
	}
}

uint64_t SPM_CRSVI_MT_NAME(_size)(void *spm)
{
	unsigned long ret;
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm_mt->spm_threads->spm;
	SPM_CRSVI_TYPE *crsvi = crsvi_mt->crsvi;

	ret = crsvi->nv*sizeof(ELEM_TYPE);
	ret += crsvi->nz*(sizeof(SPM_CRSVI_VI_TYPE) + sizeof(SPM_CRSVI_CI_TYPE));
	ret += crsvi->nrows*sizeof(SPM_CRSVI_CI_TYPE);

	return ret;
}

void SPM_CRSVI_MT_NAME(_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm_thread->spm;
	SPM_CRSVI_NAME(_destroy)(crsvi_mt->crsvi);
	free(crsvi_mt);
	free(spm_thread);
	free(spm_mt);
}

XSPMV_MT_METH_INIT(
	SPM_CRSVI_MT_NAME(_multiply),
	SPM_CRSVI_MT_NAME(_init_mmf),
	SPM_CRSVI_MT_NAME(_size),
	SPM_CRSVI_MT_NAME(_destroy),
	sizeof(ELEM_TYPE)
)

#ifdef SPM_NUMA

#include "numa.h"

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

void *SPM_CRSVI_MT_NAME(_numa_init_mmf)(char *mmf,
                                        uint64_t *rows_nr, uint64_t *cols_nr,
                                        uint64_t *nz_nr)
{
	spm_mt_t *spm_mt;
	int i;

	spm_mt = SPM_CRSVI_MT_NAME(_init_mmf)(mmf, rows_nr, cols_nr, nz_nr);

	/* sanity check */
	if (numa_available() == -1){
		perror("numa not available");
		exit(1);
	}

	/* keep a reference to original crs_vi */
	SPM_CRSVI_TYPE *crsvi = ((SPM_CRSVI_MT_TYPE *)spm_mt->spm_threads->spm)->crsvi;
	/* this will be shared */
	ELEM_TYPE *uvals = crsvi->values;
	SPM_CRSVI_CI_TYPE *rowptr = crsvi->row_ptr;
	SPM_CRSVI_CI_TYPE *colind = crsvi->col_ind;
	SPM_CRSVI_VI_TYPE *valind = crsvi->val_ind;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t  *spm_thread = spm_mt->spm_threads + i;
		SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm_thread->spm;
		/* get numa node from cpu */
		int node = numa_node_from_cpu(spm_thread->cpu);
		/* allocate new space */
		SPM_CRSVI_TYPE *numa_crsvi = numa_alloc_onnode(sizeof(SPM_CRSVI_TYPE), node);
		if (!numa_crsvi){
			perror("numa_alloc_onnode");
			exit(1);
		}
		uint64_t nnz = crsvi_mt->nnz;
		uint64_t nrows = crsvi_mt->row_end - crsvi_mt->row_start;

		numa_crsvi->col_ind = numa_alloc_onnode(sizeof(SPM_CRSVI_CI_TYPE)*nnz, node);
		numa_crsvi->row_ptr = numa_alloc_onnode(sizeof(SPM_CRSVI_CI_TYPE)*(nrows+1), node);
		numa_crsvi->val_ind = numa_alloc_onnode(sizeof(SPM_CRSVI_VI_TYPE)*(nnz), node);
		if (!numa_crsvi->col_ind || !numa_crsvi->row_ptr || !numa_crsvi->val_ind){
			perror("numa_alloc_onnode");
			exit(1);
		}
		/* copy data */
		memcpy(numa_crsvi->val_ind, valind, sizeof(SPM_CRSVI_VI_TYPE)*nnz);
		valind += nnz;
		memcpy(numa_crsvi->col_ind, colind, sizeof(SPM_CRSVI_CI_TYPE)*nnz);
		colind += nnz;
		memcpy(numa_crsvi->row_ptr, rowptr, sizeof(SPM_CRSVI_CI_TYPE)*(nrows+1));
		rowptr += nrows;
		/* make the swap */
		numa_crsvi->values = uvals;
		numa_crsvi->nrows = nrows;
		numa_crsvi->nz = nnz;
		numa_crsvi->ncols = crsvi->ncols;
		numa_crsvi->nv = crsvi->nv;
		crsvi_mt->crsvi = numa_crsvi;
	}

	/* free crsvi (except values) */
	free(crsvi->val_ind);
	free(crsvi->col_ind);
	free(crsvi->row_ptr);
	free(crsvi);

	return spm_mt;
}

void SPM_CRSVI_MT_NAME(_numa_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	int i;

	SPM_CRSVI_MT_TYPE *crsvi_mt0 = (SPM_CRSVI_MT_TYPE *)spm_mt->spm_threads->spm;
	ELEM_TYPE *uvals = crsvi_mt0->crsvi->values;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm_thread->spm;
		SPM_CRSVI_TYPE *numa_crsvi = crsvi_mt->crsvi;
		uint64_t nnz = numa_crsvi->nz;
		uint64_t nrows = crsvi_mt->row_end - crsvi_mt->row_start;
		numa_free(numa_crsvi->val_ind, sizeof(SPM_CRSVI_VI_TYPE)*nnz);
		numa_free(numa_crsvi->col_ind, sizeof(SPM_CRSVI_CI_TYPE)*nnz);
		numa_free(numa_crsvi->row_ptr, sizeof(SPM_CRSVI_CI_TYPE)*(nrows+1));
		numa_free(numa_crsvi, sizeof(SPM_CRSVI_TYPE));
	}
	free(uvals);
	free(crsvi_mt0);
	free(spm_mt->spm_threads);
	free(spm_mt);
}

uint64_t SPM_CRSVI_MT_NAME(_numa_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	uint64_t ret;
	int i;

	uint64_t uvals = ((SPM_CRSVI_MT_TYPE *)spm_mt->spm_threads->spm)->crsvi->nv;
	ret = uvals*sizeof(ELEM_TYPE);
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm_thread->spm;
		SPM_CRSVI_TYPE *numa_crsvi = crsvi_mt->crsvi;
		ret += numa_crsvi->nz*(sizeof(SPM_CRSVI_VI_TYPE) + sizeof(SPM_CRSVI_CI_TYPE));
		ret += numa_crsvi->nrows*(sizeof(SPM_CRSVI_CI_TYPE));
	}

	return ret;
}


void SPM_CRSVI_MT_NAME(_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRSVI_MT_TYPE *crsvi_mt = (SPM_CRSVI_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *values = crsvi_mt->crsvi->values;

	SPM_CRSVI_CI_TYPE *row_ptr = crsvi_mt->crsvi->row_ptr;
	SPM_CRSVI_CI_TYPE *col_ind = crsvi_mt->crsvi->col_ind;
	SPM_CRSVI_VI_TYPE *val_ind = crsvi_mt->crsvi->val_ind;
	const uint64_t row_start = crsvi_mt->row_start;
	const uint64_t row_end = crsvi_mt->row_end;
	const uint64_t nrows = row_end - row_start;
	register ELEM_TYPE yr;
	ELEM_TYPE *y = out->elements + row_start;

	unsigned long i,j;
	const register uint64_t j0 = row_ptr[0];
	for (i=0; i<nrows; i++){
		yr = (ELEM_TYPE)0;
		for (j=row_ptr[i]-j0; j<row_ptr[i+1]-j0; j++) {
			yr += (values[*(val_ind++)] * x[col_ind[j]]);
		}
		y[i] = yr;
	}
}

XSPMV_MT_METH_INIT(
 SPM_CRSVI_MT_NAME(_numa_multiply),
 SPM_CRSVI_MT_NAME(_numa_init_mmf),
 SPM_CRSVI_MT_NAME(_numa_size),
 SPM_CRSVI_MT_NAME(_numa_destroy),
 sizeof(ELEM_TYPE)
)

#endif /* SPM_NUMA  */
