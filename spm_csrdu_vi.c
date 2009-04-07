#include <stdlib.h>

#include "spm_crs_vi.h"
#include "spm_csrdu.h"
#include "spm_csrdu_vi.h"

void SPM_CSRDU_VI_NAME(_destroy)(void *m)
{
	SPM_CSRDU_VI_TYPE *csrdu_vi = (SPM_CSRDU_VI_TYPE *)m;
	free(csrdu_vi->ctl);
	free(csrdu_vi->val_ind);
	free(csrdu_vi->values);
	free(csrdu_vi);
}

uint64_t SPM_CSRDU_VI_NAME(_size)(void *m)
{
	uint64_t ret;
	SPM_CSRDU_VI_TYPE *csrdu_vi = (SPM_CSRDU_VI_TYPE *)m;
	ret = csrdu_vi->ctl_size;
	ret += csrdu_vi->nnz*(sizeof(SPM_CSRDU_VI_TYPE));
	ret += csrdu_vi->nv*(sizeof(ELEM_TYPE));

	return ret;
}

void *SPM_CSRDU_VI_NAME(_init_mmf)(char *mmf_file,
                                   uint64_t *nrows, uint64_t *ncols, uint64_t *nnz)
{
	SPM_CSRDU_TYPE *csrdu;
	SPM_CSRDU_VI_TYPE *csrdu_vi;

	csrdu = SPM_CSRDU_NAME(_init_mmf)(mmf_file, nrows, ncols, nnz);
	csrdu_vi = malloc(sizeof(SPM_CSRDU_VI_TYPE));
	if (!csrdu_vi){
		perror("malloc");
		exit(1);
	}

	csrdu_vi->nnz = csrdu->nnz;
	csrdu_vi->nrows = csrdu->nrows;
	csrdu_vi->ncols = csrdu->ncols;
	csrdu_vi->ctl = csrdu->ctl;
	csrdu_vi->ctl_size = csrdu->ctl_size;

	SPM_VI_NAME(_vals_to_vis)(
		csrdu->values, csrdu->nnz,
		&csrdu_vi->val_ind, &csrdu_vi->values, &csrdu_vi->nv);

	free(csrdu->values);
	free(csrdu);
	return csrdu_vi;
}

/* multithreading versions */

void SPM_CSRDU_VI_NAME(_mt_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt = spm_thread->spm;
	SPM_CSRDU_VI_NAME(_destroy)(csrdu_vi_mt->csrdu_vi);
	free(csrdu_vi_mt);
	free(spm_thread);
	free(spm_mt);
}

uint64_t SPM_CSRDU_VI_NAME(_mt_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt;

	csrdu_vi_mt = (SPM_CSRDU_VI_MT_TYPE *)spm_mt->spm_threads->spm;
	return SPM_CSRDU_VI_NAME(_size)(csrdu_vi_mt->csrdu_vi);
}

void *SPM_CSRDU_VI_NAME(_mt_init_mmf)(char *mmf_file,
                                      uint64_t *nrows, uint64_t *ncols, uint64_t *nnz)
{
	SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt;
	SPM_CSRDU_VI_TYPE *csrdu_vi;
	SPM_CSRDU_TYPE *csrdu;
	spm_csrdu_mt_t *csrdu_mt;
	spm_mt_t *spm_mt;
	int nr_threads, i;

	spm_mt = SPM_CSRDU_NAME(_mt_init_mmf)(mmf_file, nrows, ncols, nnz);
	nr_threads = spm_mt->nr_threads;

	csrdu_vi = malloc(sizeof(SPM_CSRDU_VI_TYPE));
	csrdu_vi_mt = malloc(sizeof(SPM_CSRDU_VI_MT_TYPE)*nr_threads);
	if (!csrdu_vi_mt || !csrdu_vi){
		perror("malloc");
		exit(1);
	}

	csrdu_mt = spm_mt->spm_threads->spm;
	csrdu = csrdu_mt->csrdu;

	/* csrdu to csrdu_vi */
	csrdu_vi->nnz = csrdu->nnz;
	csrdu_vi->nrows = csrdu->nrows;
	csrdu_vi->ncols = csrdu->ncols;
	csrdu_vi->ctl = csrdu->ctl;
	csrdu_vi->ctl_size = csrdu->ctl_size;
	SPM_VI_NAME(_vals_to_vis)(
		csrdu->values, csrdu->nnz,
		&csrdu_vi->val_ind, &csrdu_vi->values, &csrdu_vi->nv);

	for (i=0; i < nr_threads; i++){
		csrdu_vi_mt[i].nnz = csrdu_mt[i].nnz;
		csrdu_vi_mt[i].row_start = csrdu_mt[i].row_start;
		csrdu_vi_mt[i].val_start = csrdu_mt[i].val_start;
		csrdu_vi_mt[i].ctl_start = csrdu_mt[i].ctl_start;
		csrdu_vi_mt[i].csrdu_vi = csrdu_vi;
		spm_mt->spm_threads[i].spm = csrdu_vi_mt + i;
	}

	free(csrdu_mt);
	free(csrdu->values);
	free(csrdu);

	return spm_mt;
}

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

static void check_numa()
{
	if (numa_available() == -1){
		perror("numa not available");
		exit(1);
	}
}

void *SPM_CSRDU_VI_NAME(_mt_numa_init_mmf)(char *mmf_file,
                                           uint64_t *nrows_ptr, uint64_t *ncols_ptr,
                                           uint64_t *nnz_ptr)
{
	spm_mt_t *spm_mt;
	int i;
	SPM_CSRDU_VI_TYPE *csrdu_vi;

	check_numa();

	spm_mt = SPM_CSRDU_VI_NAME(_mt_init_mmf)(mmf_file, nrows_ptr, ncols_ptr, nnz_ptr);
	/* keep a reference to the original csrdu_vi */
	csrdu_vi = ((SPM_CSRDU_VI_MT_TYPE *)spm_mt->spm_threads[0].spm)->csrdu_vi;
	ELEM_TYPE *values = csrdu_vi->values;
	SPM_CRSVI_VI_TYPE *val_ind = csrdu_vi->val_ind;
	uint8_t *ctl = csrdu_vi->ctl;
	int nr_threads = spm_mt->nr_threads;

	for (i=0; i<nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt = (SPM_CSRDU_VI_MT_TYPE *)spm_thread->spm;
		/* get numa node from cpu */
		int node = numa_node_from_cpu(spm_thread->cpu);
		/* allocate new space */
		SPM_CSRDU_VI_TYPE *numa_csrdu_vi = numa_alloc_onnode(sizeof(SPM_CSRDU_VI_TYPE), node);
		if (!numa_csrdu_vi){
			perror("numa_alloc_onnode");
			exit(1);
		}
		uint64_t ctl_start = csrdu_vi_mt->ctl_start;
		uint64_t ctl_end = (i < nr_threads -1) ? (csrdu_vi_mt+1)->ctl_start : csrdu_vi->ctl_size;
		uint64_t ctl_size = ctl_end - ctl_start;
		uint64_t row_start = csrdu_vi_mt->row_start;
		uint64_t row_end = (i < nr_threads -1) ? (csrdu_vi_mt+1)->row_start : csrdu_vi->nrows;
		uint64_t nnz = csrdu_vi_mt->nnz;
		numa_csrdu_vi->val_ind = numa_alloc_onnode(sizeof(SPM_CRSVI_VI_TYPE)*nnz, node);
		/*
		 * Note that the newly allocated ctl should have the same alignment with
		 * the previous part
		 */
		unsigned long align = (unsigned long)(ctl + ctl_start) & (8UL-1UL);
		numa_csrdu_vi->ctl = numa_alloc_onnode(ctl_size + align, node);
		if (!numa_csrdu_vi->val_ind || !numa_csrdu_vi->ctl){
			perror("numa_alloc_onnode");
			exit(1);
		}
		/* copy data */
		memcpy(numa_csrdu_vi->val_ind, val_ind, sizeof(SPM_CRSVI_VI_TYPE)*nnz);
		val_ind += nnz;
		memcpy(numa_csrdu_vi->ctl + align, ctl + ctl_start, ctl_size);
		/* make the swap */
		numa_csrdu_vi->nnz = nnz;
		numa_csrdu_vi->ncols = csrdu_vi->ncols;
		numa_csrdu_vi->nv = csrdu_vi->nv;
		numa_csrdu_vi->nrows = row_end - row_start;
		numa_csrdu_vi->ctl_size = ctl_size;
		numa_csrdu_vi->values = values;
		csrdu_vi_mt->csrdu_vi = numa_csrdu_vi;
		csrdu_vi_mt->ctl_start = align;
	}

	free(csrdu_vi->ctl);
	free(csrdu_vi->val_ind);
	free(csrdu_vi);

	return spm_mt;
}

void SPM_CSRDU_VI_NAME(_mt_numa_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	int i;
	SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt0 = (SPM_CSRDU_VI_MT_TYPE *)spm_mt->spm_threads->spm;
	ELEM_TYPE *uvals = csrdu_vi_mt0->csrdu_vi->values;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt = (SPM_CSRDU_VI_MT_TYPE *)spm_thread->spm;
		SPM_CSRDU_VI_TYPE *csrdu_vi = csrdu_vi_mt->csrdu_vi;
		numa_free(csrdu_vi->val_ind, csrdu_vi->nnz);
		numa_free(csrdu_vi->ctl, csrdu_vi->ctl_size);
		numa_free(csrdu_vi, sizeof(SPM_CSRDU_VI_TYPE));
	}
	free(uvals);
	free(csrdu_vi_mt0);
	free(spm_mt->spm_threads);
	free(spm_mt);
}

uint64_t SPM_CSRDU_VI_NAME(_mt_numa_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	uint64_t ret;
	int i;

	uint64_t uvals = ((SPM_CSRDU_VI_MT_TYPE *)spm_mt->spm_threads->spm)->csrdu_vi->nv;
	ret = uvals*sizeof(ELEM_TYPE);
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CSRDU_VI_MT_TYPE *csrdu_vi_mt = (SPM_CSRDU_VI_MT_TYPE *)spm_thread->spm;
		SPM_CSRDU_VI_TYPE *numa_csrdu_vi = csrdu_vi_mt->csrdu_vi;
		ret += numa_csrdu_vi->nnz*(sizeof(SPM_CRSVI_VI_TYPE));
		ret += numa_csrdu_vi->ctl_size;
	}

	return ret;
}

#endif /* SPM_NUMA */
