#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <inttypes.h>

#ifndef SPM_CRSVH_CI_BITS
#define SPM_CRSVH_CI_BITS 32
#endif

#define SPM_CRSVH_CI_TYPE UINT_TYPE(SPM_CRSVH_CI_BITS)
#define SPM_CRS_BITS SPM_CRSVH_CI_BITS

#include "phash.h"
#include "vector.h"
#include "huffman.h"
#include "bitutils.h"
#include "dynarray.h"
#include "mmf.h"
#include "spmv_method.h"

#include "spm_mt.h"
#include "spm_crs.h"
#include "spm_crs_mt.h"
#include "spm_crs_vh_mt.h"

#define CRSVH_ENV_LIMIT "CRSVH_LIMIT"
#define CRSVH_ENV_BITS "CRSVH_BITS"

static void crsvh_mt_compress_vals(spm_mt_thread_t *spm_thread, SPM_CRS_VH_MT_TYPE *crsvh,
                                   int limit, int b)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
	SPM_CRS_TYPE *crs = crs_mt->crs;

	SPM_CRSVH_CI_TYPE *row_ptr = crs->row_ptr;
	unsigned long row_s = crs_mt->row_start;
	unsigned long row_e = crs_mt->row_end;
	void *vals = (void *)(crs->values + row_ptr[row_s]);
	unsigned long vals_s = (row_ptr[row_e] -row_ptr[row_s])*sizeof(ELEM_TYPE);

	phash_t *syms;
	int bits=0;
	huff_mksymcodes_split(vals, vals_s, limit, b, &syms, &crsvh->htree, &bits);
	//printf("crsvh_mt: encoding: limit:%d bits:%d rbits:%d\n", limit, b, bits);
	do_huff_encode(vals, vals_s, &crsvh->hvals, &crsvh->hvals_bits, syms, bits);
	phash_free(syms);

	crsvh->nz = vals_s/sizeof(ELEM_TYPE);
	crsvh->row_start = crs_mt->row_start;
	crsvh->row_end = crs_mt->row_end;
	crsvh->nrows = crs_mt->crs->nrows;
	crsvh->ncols = crs->ncols;
	crsvh->col_ind = crs->col_ind;
	crsvh->row_ptr = crs->row_ptr;
	printf("cpu:%d nz:%lu row_start:%lu row_end:%lu\n",
	        spm_thread->cpu,
		crsvh->nz,
		crsvh->row_start,
		crsvh->row_end);

	spm_thread->spm = crsvh;
}

spm_mt_t *SPM_CRS_VH_MT_NAME(_init_mmf)(char *mmf_file,
                                        unsigned long *rows_nr, unsigned long *cols_nr,
                                        unsigned long *nz_nr)
{
	spm_mt_t *spm_mt;

	void *s;
	int limit=0, bits=64;
	if ( (s = getenv(CRSVH_ENV_LIMIT)) ){
		limit = atol(s);
	}
	if ( (s = getenv(CRSVH_ENV_BITS)) ){
		bits = atol(s);
		assert((bits==8) || (bits==16) || (bits==32) || (bits==64));
	}

	spm_mt = SPM_CRS_MT_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr);

	SPM_CRS_VH_MT_TYPE *crsvh;
	crsvh = malloc(sizeof(SPM_CRS_VH_MT_TYPE)*spm_mt->nr_threads);
	if (!crsvh){
		perror("malloc");
		exit(1);
	}

	int i;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_mt->spm_threads[0].spm;
	for (i=0; i < spm_mt->nr_threads; i++){
		spm_mt_thread_t *spm_thr = &spm_mt->spm_threads[i];
		crsvh_mt_compress_vals(spm_thr, &crsvh[i], limit, bits);
	}
	free(crs_mt[0].crs->values);
	free(crs_mt[0].crs);
	free(crs_mt);

	return spm_mt;
}

unsigned long SPM_CRS_VH_MT_NAME(_size)(spm_mt_t *spm_mt)
{
	unsigned long ret=0, nz=0, nrows = 0;
	int i;
	for (i=0; i < spm_mt->nr_threads; i++){
		SPM_CRS_VH_MT_TYPE *crsvh = spm_mt->spm_threads[i].spm;
		ret += crsvh->hvals_bits/8 + (crsvh->hvals_bits % 8);
		nz += crsvh->nz;
		nrows = crsvh->nrows;
	}
	ret += nz*(sizeof(SPM_CRSVH_CI_TYPE));
	ret += nrows*sizeof(SPM_CRSVH_CI_TYPE);
	return ret;
}

void SPM_CRS_VH_MT_NAME(_multiply) (void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_VH_MT_TYPE *crs_vh = (SPM_CRS_VH_MT_TYPE *)spm;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *x = in->elements;
	const register SPM_CRSVH_CI_TYPE *row_ptr = crs_vh->row_ptr;
	const register SPM_CRSVH_CI_TYPE *col_ind = crs_vh->col_ind;
	const register unsigned long rs = crs_vh->row_start;
	const register unsigned long re = crs_vh->row_end;
	register huff_node_t *tree = crs_vh->htree, *n;
	register unsigned char *hvals = crs_vh->hvals;
	register unsigned long bitcnt = 0;
	register ELEM_TYPE yr;
	register ELEM_TYPE val;

	register unsigned long i, j=0;
	for(i=rs; i<re; i++) {
		yr = (ELEM_TYPE)0;
		__asm__ __volatile__ ("# loop start\n\t");
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) {
			/* get value */
			for (n=tree;;){
				n = bit_test(hvals, bitcnt) ? n->r : n->l;
				bitcnt++;
				if (n->flags & HUFF_NODE_LEAF_FL){
					val = *((ELEM_TYPE *)((&n->val))); /* XXX: This won't work for floats */
					break;
				}
			}
			yr += (val * x[col_ind[j]]);
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
	SPM_CRS_VH_MT_NAME(_multiply),
	SPM_CRS_VH_MT_NAME(_init_mmf),
	SPM_CRS_VH_MT_NAME(_size),
	sizeof(ELEM_TYPE)
)
