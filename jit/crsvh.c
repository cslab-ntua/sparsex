#include <stdlib.h>

#include "method.h"
#include "huff_llvm_jit.h"

#define SPM_CRSVH_CI_BITS CI_BITS

#include "spm_mt.h"
#include "spm_crs_vh.h"
#include "spm_crs_vh_mt.h"
#include "jit/crsvh.h"
#include "jit/vh_jit.h"

#define VHJIT_TYPE CON3(spm_crs, CI_BITS, _vhjit_t)
#define _NAME(name) CON5(spm_crs, CI_BITS, _vhjit_, ELEM_TYPE, name)

#if ELEM_TYPE == double
	#define OUT_SYMBOL   HJ_OUT_DOUBLE
	#define JIT_FN_TYPE  huff_decoder_double_fn_t
#elif ELEM_TYPE == float
	#define OUT_SYMBOL   HJ_OUT_FLOAT
	#define JIT_FN_TYPE  huff_decoder_float_fn_t
#endif

method_t *_NAME(_init_method)(char *mmf_file,
                              unsigned long *rows_nr, unsigned long *cols_nr,
                              unsigned long *nz_nr)
{
	method_t *method = malloc(sizeof(*method));
	if (!method){
		perror("malloc");
		exit(1);
	}

	SPM_CRSVH_TYPE *crsvh = SPM_CRS_VH_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);
	VHJIT_TYPE *vhjit = malloc(sizeof(*vhjit));
	if (!vhjit){
		perror("malloc");
		exit(1);
	}

	vhjit->col_ind = crsvh->col_ind;
	vhjit->row_ptr = crsvh->row_ptr;
	vhjit->nz = crsvh->nz;
	vhjit->nrows = crsvh->nrows;
	vhjit->ncols = crsvh->ncols;
	vhjit->hvals = crsvh->hvals;
	vhjit->row_start = 0;
	vhjit->row_end = crsvh->nrows;

	printf("csrvh symbol bits: %lu\n", crsvh->hsym_bits);
	void *hj = vhjit_init(crsvh->htree, crsvh->hsym_bits, OUT_SYMBOL);

	method->name = "spmv_vhjit_blah";
	method->fn = vhjit_get_spmvfn(hj, CI_BITS);
	method->next = NULL;
	method->data = vhjit;

	free(crsvh);
	return method;
}

#define _NAME_MT(name) CON5(spm_crs, CI_BITS, _vhjit_mt_, ELEM_TYPE, name)

spm_mt_t *_NAME_MT(_init_method)(char *mmf_file,
                                 unsigned long *rows_nr, unsigned long *cols_nr,
                                 unsigned long *nz_nr)
{
	spm_mt_t *spm_mt;
	spm_mt = SPM_CRS_VH_MT_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);

	VHJIT_TYPE *vhjit_mt = malloc(sizeof(*vhjit_mt)*spm_mt->nr_threads);
	if (!vhjit_mt){
		perror("malloc");
		exit(1);
	}

	void *hj[spm_mt->nr_threads];

	int i;
	SPM_CRS_VH_MT_TYPE *crsvh_mt = spm_mt->spm_threads[0].spm;
	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thr = &spm_mt->spm_threads[i];
		SPM_CRS_VH_MT_TYPE *crsvh_mt = (SPM_CRS_VH_MT_TYPE *)spm_thr->spm;

		vhjit_mt[i].row_ptr = crsvh_mt->row_ptr;
		vhjit_mt[i].col_ind = crsvh_mt->col_ind;
		vhjit_mt[i].nz = crsvh_mt->nz;
		vhjit_mt[i].ncols = crsvh_mt->ncols;
		vhjit_mt[i].nrows = crsvh_mt->nrows;
		vhjit_mt[i].row_start = crsvh_mt->row_start;
		vhjit_mt[i].row_end = crsvh_mt->row_end;
		vhjit_mt[i].hvals = crsvh_mt->hvals;
		vhjit_mt[i].hvals_bits = crsvh_mt->hvals_bits;

		hj[i] = vhjit_init(crsvh_mt->htree, crsvh_mt->hsym_bits, OUT_SYMBOL);

		spm_thr->spm = &vhjit_mt[i];
	}
	free(crsvh_mt);

	for (i=0; i<spm_mt->nr_threads; i++){
		spm_mt->spm_threads[i].spmv_fn = vhjit_get_spmvfn(hj[i], CI_BITS);
	}

	return spm_mt;
}


unsigned long _NAME_MT(_size)(spm_mt_t *spm_mt)
{
	unsigned long ret=0, nz=0, nrows = 0;
	int i;
	for (i=0; i < spm_mt->nr_threads; i++){
		VHJIT_TYPE *vhjit = spm_mt->spm_threads[i].spm;
		ret += vhjit->hvals_bits/8 + (vhjit->hvals_bits % 8);
		nz += vhjit->nz;
		nrows = vhjit->nrows;
	}
	ret += nz*(CI_BITS/8);
	ret += nrows*(CI_BITS/8);
	return ret;
}
