#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <inttypes.h>

#include "macros.h"

#ifndef SPM_CRSVH_CI_BITS
#define SPM_CRSVH_CI_BITS 32
#endif

#define SPM_CRSVH_CI_TYPE UINT_TYPE(SPM_CRSVH_CI_BITS)

#include "phash.h"
#include "vector.h"
#include "huffman.h"
#include "bitutils.h"
#include "dynarray.h"
#include "mmf.h"
#include "spmv_method.h"

#define SPM_CRS_BITS SPM_CRSVH_CI_BITS

#include "spm_crs.h"
#include "spm_crs_vh.h"

#define CRSVH_ENV_LIMIT "CRSVH_LIMIT"
#define CRSVH_ENV_BITS "CRSVH_BITS"

static void crsvh_compress_vals(SPM_CRSVH_TYPE *crs_vh, SPM_CRS_TYPE *crs)
{
	phash_t *syms;
	int limit=0, b=64, bits=0;
	unsigned long vals_s = crs->nz*sizeof(ELEM_TYPE);

	void *s;
	if ( (s = getenv(CRSVH_ENV_LIMIT)) ){
		limit = atol(s);
	}
	if ( (s = getenv(CRSVH_ENV_BITS)) ){
		b = atol(s);
		assert((b==8) || (b==16) || (b==32) || (b==64));
	}

	huff_mksymcodes_split((void *)crs->values, vals_s, limit, b, &syms, &crs_vh->htree, &bits);
	printf("crsvh: encoding: limit:%d bits:%d rbits:%d\n", limit, b, bits);
	do_huff_encode(crs->values, vals_s,
	               &crs_vh->hvals, &crs_vh->hvals_bits,
		       syms, bits);

	phash_free(syms);
	free(crs->values);
	free(crs);
}


SPM_CRSVH_TYPE *SPM_CRS_VH_NAME(_init_mmf) (char *mmf_file,
                                       unsigned long *rows_nr, unsigned long *cols_nr,
                                       unsigned long *nz_nr)
{
	SPM_CRSVH_TYPE *crs_vh;
	SPM_CRS_TYPE *crs;

	crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);
	crs_vh = malloc(sizeof(SPM_CRSVH_TYPE));
	if (!crs_vh){
		perror("malloc");
		exit(1);
	}

	crs_vh->col_ind = crs->col_ind;
	crs_vh->row_ptr = crs->row_ptr;
	crs_vh->nz = crs->nz;
	crs_vh->nrows = crs->nrows;
	crs_vh->ncols = crs->ncols;

	crsvh_compress_vals(crs_vh, crs);

	return crs_vh;
}

void SPM_CRS_VH_NAME(_destroy)(SPM_CRSVH_TYPE *crs_vh)
{
	free(crs_vh->hvals);
	free(crs_vh->col_ind);
	free(crs_vh->row_ptr);
	free(crs_vh);
}

unsigned long SPM_CRS_VH_NAME(_size)(SPM_CRSVH_TYPE *crsvh)
{
	unsigned long ret;
	ret  = crsvh->hvals_bits/8 + 1*(crsvh->hvals_bits % 8);
	ret += crsvh->nz*(sizeof(SPM_CRSVH_CI_TYPE));
	ret += crsvh->ncols*sizeof(SPM_CRSVH_CI_TYPE);

	return ret;
}

void SPM_CRS_VH_NAME(_multiply) (void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRSVH_TYPE *crs_vh = (SPM_CRSVH_TYPE *)spm;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *x = in->elements;
	const register SPM_CRSVH_CI_TYPE *row_ptr = crs_vh->row_ptr;
	const register SPM_CRSVH_CI_TYPE *col_ind = crs_vh->col_ind;
	register huff_node_t *tree = crs_vh->htree, *n;
	register unsigned char *hvals = crs_vh->hvals;
	const register unsigned long nrows = crs_vh->nrows;
	register unsigned long bitcnt = 0;
	register ELEM_TYPE yr;
	register ELEM_TYPE val;

	register unsigned long i, j=0;
	for(i=0; i<nrows; i++) {
		yr = (ELEM_TYPE)0;
		__asm__ __volatile__ ("# loop start\n\t");
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) {
			/* get value */
			for (n=tree;;){
				n = bit_test(hvals, bitcnt) ? n->r : n->l;
				bitcnt++;
				if (n->flags & HUFF_NODE_LEAF_FL){
					val = *((ELEM_TYPE *)(&n->val)); /* XXX: This won't work for floats */
					break;
				}
			}
			yr += (val * x[col_ind[j]]);
			//printf("i:%lu j:%lu val:%lf col_ind:%lu x:%lf yr:%lf\n", i, j, val, (unsigned long)col_ind[j], x[col_ind[j]], yr);
		}
		y[i] = yr;
		__asm__ __volatile__ ("# loop end\n\t");
	}
	#if 0
	if (bitcnt != crs_vh->hvals_bits){
		printf("bitcnt: %lu hvals_bits: %lu\n", bitcnt, crs_vh->hvals_bits);
		assert(0);
	}
	#endif
}

#define XSPMV_METH_INIT(x,y,z,w) SPMV_METH_INIT(x,y,z,w)
XSPMV_METH_INIT(
	SPM_CRS_VH_NAME(_multiply),
	SPM_CRS_VH_NAME(_init_mmf),
	SPM_CRS_VH_NAME(_size),
	sizeof(ELEM_TYPE)
)
