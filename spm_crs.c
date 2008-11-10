#include <stdlib.h>
#include <emmintrin.h>
#include <inttypes.h>
#include <assert.h>

#include "macros.h"
	
#ifndef SPM_CRS_BITS
#define SPM_CRS_BITS 64
#endif

#define SPM_CRS_IDX_TYPE UINT_TYPE(SPM_CRS_BITS)

#include "vector.h"
#include "spm_crs.h"
#include "dynarray.h"
#include "mmf.h"
#include "method.h"
#include "spm_parse.h"

#define SPM_CRS_NAME(name) CON5(spm_crs, SPM_CRS_BITS, _, ELEM_TYPE, name)
#define SPM_CRS_TYPE SPM_CRS_NAME(_t)

typedef SPM_CRS_IDX_TYPE spm_crs_idx_t;
typedef SPM_CRS_TYPE spm_crs_t;

struct crs_state_s {
	spm_crs_t *crs;
	dynarray_t *sp_values, *sp_colind, *sp_rowptr;
};
typedef struct crs_state_s crs_state_t;

static void _initialize(parse_state_t *parse_state, void **crs_state_ptr,
                        unsigned long rows_nr, unsigned long cols_nr,
                        unsigned long nz_nr)
{
	spm_crs_t *crs;
	crs_state_t *crs_state;

	crs = malloc(sizeof(spm_crs_t));
	crs_state = malloc(sizeof(crs_state_t));
	if ( !crs || !crs_state ){
		perror("spm_crs_init: malloc");
		exit(1);
	}

	crs->ncols = cols_nr;
	crs->nrows = rows_nr;
	crs->nz = 0;

	crs_state->crs = crs;
	crs_state->sp_values = dynarray_create(sizeof(ELEM_TYPE), nz_nr);
	crs_state->sp_colind = dynarray_create(sizeof(spm_crs_idx_t), nz_nr);
	crs_state->sp_rowptr = dynarray_create(sizeof(spm_crs_idx_t), rows_nr);
	spm_crs_idx_t *rowptr = dynarray_alloc(crs_state->sp_rowptr);
	*rowptr = 0;
	assert(parse_state->row == 0);
	*crs_state_ptr = crs_state;
}

static void *_finalize(parse_state_t *parse_state, void *crs_st)
{
	crs_state_t *crs_state = crs_st;
	spm_crs_t *crs = crs_state->crs;
	spm_crs_idx_t *rowptr = dynarray_alloc(crs_state->sp_rowptr);

	assert(crs->nz == dynarray_size(crs_state->sp_values));
	*rowptr = crs->nz;
	crs->values = dynarray_destroy(crs_state->sp_values);
	crs->col_ind = dynarray_destroy(crs_state->sp_colind);
	crs->row_ptr = dynarray_destroy(crs_state->sp_rowptr);

	free(crs_state);

	return crs;
}

static void _finalize_unit(parse_state_t *parse_state, void *crs_st)
{
	crs_state_t *crs_state = crs_st;
	unsigned long l;
	ELEM_TYPE *vals = dynarray_alloc_nr(crs_state->sp_values, parse_state->cur_values_idx);
	spm_crs_idx_t *colind = dynarray_alloc_nr(crs_state->sp_colind, parse_state->cur_values_idx);
	for ( l=0; l < parse_state->cur_values_idx; l++){
		vals[l] = (ELEM_TYPE)parse_state->cur_values[l];
		colind[l] = parse_state->cur_unit_col + l;
	}
	crs_state->crs->nz += parse_state->cur_values_idx;
}

static void _finalize_row(parse_state_t *parse_state, void *crs_st)
{
	crs_state_t *crs_state = crs_st;
	unsigned long bsize = dynarray_size(crs_state->sp_rowptr);
	unsigned long empty_rows = parse_state->cur_row - bsize + 1;
	spm_crs_idx_t *rowptr = dynarray_alloc_nr(crs_state->sp_rowptr, 1+empty_rows);
	unsigned long i;
	for ( i=0; i < empty_rows; i++){
		*(rowptr+i) = *(rowptr+i-1);
	}
	rowptr[i] = crs_state->crs->nz;
}

static spm_parser_t SPM_CRS_NAME(_parser) = {
	.initialize = _initialize,
	.finalize = _finalize,
	.finalize_unit = _finalize_unit,
	.finalize_row = _finalize_row
};

SPM_CRS_TYPE *SPM_CRS_NAME(_init_mmf) (char *mmf_file, 
                                       unsigned long *rows_nr, unsigned long *cols_nr,
                                       unsigned long *nz_nr)
{
	SPM_CRS_TYPE *crs;
	FILE *f;
	f = mmf_init(mmf_file, rows_nr, cols_nr, nz_nr);
	crs  = spm_parse_seq(f, *rows_nr, *cols_nr, *nz_nr, &SPM_CRS_NAME(_parser));
	fclose(f);

	return crs;
}

void SPM_CRS_NAME(_destroy)(SPM_CRS_TYPE *crs)
{
	free(crs->values);
	free(crs->col_ind);
	free(crs->row_ptr);
	free(crs);
}

void SPM_CRS_NAME(_multiply) (void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	//printf("++%s\n", __FUNCTION__);
	SPM_CRS_TYPE *crs = (SPM_CRS_TYPE *)spm;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *values = crs->values;
	SPM_CRS_IDX_TYPE *row_ptr = crs->row_ptr;
	SPM_CRS_IDX_TYPE *col_ind = crs->col_ind;
	unsigned long n = crs->nrows;
	register ELEM_TYPE yr;

	unsigned long i, j;
	for(i=0; i<n; i++) {
		yr = (ELEM_TYPE)0;
		//printf("row_ptr_i: %lu row_ptr_i+1: %lu \n", (unsigned long)row_ptr[i], (unsigned long)row_ptr[i+1]);
		for(j=row_ptr[i]; j<row_ptr[i+1]; j++) { 
			yr += (values[j] * x[col_ind[j]]);
		}
		y[i] = yr;
		//printf("++y[%lu] = %lf\n", i, yr);
		#if 0
		__asm__ __volatile__ (
			" movntq %[val], (%[mem]) \n\t" 
			: 
			: [val] "x" (yr), [mem] "r" (y+i)
		);
		#endif
	}
}

/*
#define XMETHOD_INIT(x,y) METHOD_INIT(x,y)
XMETHOD_INIT(SPM_CRS_NAME(_multiply), SPM_CRS_NAME(_init_mmf))
*/

#if 0
#include "spmv_method.h"
#define XSPMV_METH_INIT(x,y,z) SPMV_METH_INIT(x,y,z)
XSPMV_METH_INIT(
	SPM_CRS_NAME(_multiply), 
	SPM_CRS_NAME(_init_mmf), 
	SPM_CRS_NAME(_papaki)
)
#endif
