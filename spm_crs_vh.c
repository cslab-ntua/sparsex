#include <stdlib.h>
#include <assert.h>
#include <time.h>

#if 0
#include <inttypes.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#endif

#ifndef SPM_CRSVH_CI_BITS
#define SPM_CRSVH_CI_BITS 32
#endif

#if SPM_CRSVH_CI_BITS == 32
#define  SPM_CRSVH_CI_TYPE uint32_t
#elif SPM_CRSVH_CI_BITS == 64
#define SPM_CRSVH_CI_TYPE uint64_t
#else
#error "SPM_CRSHVH_CI_BITS not 32 or 64"
#endif

#if 0
#include "vector.h"
#include "ext_prog.h"
#include "spmv_method.h"
#include "phash.h"
#endif

#include "phash.h"
#include "huffman.h"
#include "bitutils.h"
#include "dynarray.h"
#include "mmf.h"

#include "spm_crs_vh.h"

#define _CON5(a,b,c,d,e) a ## b ## c ## d ## e
#define CON5(a,b,c,d,e) _CON5(a,b,c,d,e)
#define SPM_CRS_VH_NAME(name) \
	CON5(spm_crs, SPM_CRSVH_CI_BITS, _vh_, ELEM_TYPE, name)
#define SPM_CRSVH_TYPE SPM_CRS_VH_NAME(_t)

struct crs_vh_state_s {
	SPM_CRSVH_TYPE *crs_vh;
	dynarray_t     *sp_colind, *sp_rowptr, *sp_values;
};
typedef struct crs_vh_state_s crs_vh_state_t;

static void crsvh_initialize(crs_vh_state_t **crs_vh_state_ptr,
                             unsigned long rows_nr, unsigned long cols_nr,
                             unsigned long nz_nr)
{
	SPM_CRSVH_TYPE *crsvh;
	crs_vh_state_t *crsvh_st;

	crsvh = malloc(sizeof(SPM_CRSVH_TYPE));
	crsvh_st = malloc(sizeof(crs_vh_state_t));
	if ( !crsvh || !crsvh_st ){
		perror("spm_crs_vh_init: malloc");
		exit(1);
	}

	crsvh->ncols = cols_nr;
	crsvh->nrows = rows_nr;
	crsvh->nz = 0;

	crsvh_st->crs_vh = crsvh;

	crsvh_st->sp_values = dynarray_create(sizeof(ELEM_TYPE), 4096);
	crsvh_st->sp_colind = dynarray_create(sizeof(SPM_CRSVH_CI_TYPE), nz_nr);
	crsvh_st->sp_rowptr = dynarray_create(sizeof(SPM_CRSVH_CI_TYPE), rows_nr);

	/* if this fails we need to make phash code more generic */
	assert(sizeof(double) == sizeof(unsigned long));
	*crs_vh_state_ptr = crsvh_st;
}

#define CRSVH_ENV_LIMIT "CRSVH_LIMIT"
#define CRSVH_ENV_BITS "CRSVH_BITS"

static void crsvh_compress(crs_vh_state_t *crsvh_st)
{
	SPM_CRSVH_TYPE *crs_vh = crsvh_st->crs_vh;
	unsigned long vals_s = crs_vh->nz*sizeof(ELEM_TYPE);
	void *vals = dynarray_destroy(crsvh_st->sp_values);
	phash_t *syms;
	int limit=0, b=64, bits=0;
	void *s;

	if ( (s = getenv(CRSVH_ENV_LIMIT)) ){
		limit = atol(s);
	}
	if ( (s = getenv(CRSVH_ENV_BITS)) ){
		b = atol(s);
		assert((b==8) || (b==16) || (b==32) || (b==64));
	}

	huff_mksymcodes_split(vals, vals_s, limit, b, &syms, &crs_vh->htree, &bits);
	printf("crsvh: encoding: limit:%d bits:%d rbits:%d\n", limit, b, bits);
	do_huff_encode(vals, vals_s,
	               &crs_vh->hvals, &crs_vh->hvals_bits,
		       syms, bits);
	phash_free(syms);
	free(vals);
}

static void *crsvh_finalize(crs_vh_state_t *crsvh_st)
{
	SPM_CRSVH_TYPE *crs_vh = crsvh_st->crs_vh;
	SPM_CRSVH_CI_TYPE *rowptr = dynarray_alloc(crsvh_st->sp_rowptr);

	assert(crs_vh->nz == dynarray_size(crsvh_st->sp_values));
	*rowptr = crs_vh->nz;

	crsvh_compress(crsvh_st);
	crs_vh->col_ind = dynarray_destroy(crsvh_st->sp_colind);
	crs_vh->row_ptr = dynarray_destroy(crsvh_st->sp_rowptr);

	free(crsvh_st);

	return crs_vh;
}

SPM_CRSVH_TYPE *SPM_CRS_VH_NAME(_init_mmf) (char *mmf_file,
                                       unsigned long *rows_nr, unsigned long *cols_nr,
                                       unsigned long *nz_nr)
{
	crs_vh_state_t *crsvh_st;
	FILE *f;
	double val;
	unsigned long row, col, prev_row=0, i;
	time_t t0, tn;
	char *report;
	unsigned long report_rows;

	f = mmf_init(mmf_file, rows_nr, cols_nr, nz_nr);

	crsvh_initialize(&crsvh_st, *rows_nr, *cols_nr, *nz_nr);
	SPM_CRSVH_TYPE *crsvh = crsvh_st->crs_vh;

	SPM_CRSVH_CI_TYPE *rowptr = dynarray_alloc(crsvh_st->sp_rowptr);
	*rowptr = 0;
	prev_row = 0;

	report = getenv("SPM_CRSVH_REPORT");
	if (report != NULL){
		report_rows = atol(report);
		if (!report_rows)
			report_rows = 1024;
	} else {
		report_rows = 0;
	}

	unsigned long empty_rows;
	t0 = time(NULL);
	while (mmf_get_next(f, &row, &col, &val)){
		crsvh->nz++;

		/* row indices */
		if (prev_row < row){
			#if 1
			if (report_rows && row % report_rows == 0){
				tn = time(NULL);
				double ratio = (double)row/(tn - t0);
				unsigned long remaining = *rows_nr - row;
				printf("%s [ %lf m] remaining: %lu rows/sec:%lf ETA:%lf m\n",
				        mmf_file, (double)(tn-t0)/60.0, remaining, ratio, (double)remaining/(ratio*60.0));
			}
			#endif
			empty_rows = row -prev_row -1;
			rowptr = dynarray_alloc_nr(crsvh_st->sp_rowptr,1+empty_rows);
			for (i=0; i<empty_rows; i++){
				*(rowptr+i) = *(rowptr+i-1);
			}
			*(rowptr+i) = crsvh->nz - 1;
			prev_row = row;
		}

		/* column indices */
		SPM_CRSVH_CI_TYPE *colind = dynarray_alloc(crsvh_st->sp_colind);
		*colind = col;

		/* values */
		ELEM_TYPE *v = dynarray_alloc(crsvh_st->sp_values);
		*v = (ELEM_TYPE)val;
	}

	crsvh_finalize(crsvh_st);
	fclose(f);

	return crsvh;
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
