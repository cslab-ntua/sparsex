#include <stdlib.h>

#include "method.h"
#include "huff_llvm_jit.h"

#define SPM_CRSVH_CI_BITS CI_BITS

#include "spm_crs_vh.h"
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

	printf("csrvh symbol bits: %lu\n", crsvh->hsym_bits);
	void *hj = vhjit_init(crsvh->htree, crsvh->hsym_bits, OUT_SYMBOL);

	method->name = "spmv_vhjit_blah";
	method->fn = vhjit_get_spmvfn(hj, CI_BITS);
	method->next = NULL;
	method->data = vhjit;

	free(crsvh);
	return method;
}

#define VHJIT_MT_TYPE CON3(spm_crs, CI_BITS, _vhjit_mt_t)
#define _NAME_MT(name) CON5(spm_crs, CI_BITS, _vhjit_mt_, ELEM_TYPE, name)
