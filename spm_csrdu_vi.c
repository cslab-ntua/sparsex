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
