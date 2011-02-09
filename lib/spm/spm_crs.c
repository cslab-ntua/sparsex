/*
 * spm_crs.c
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdlib.h>
//#if defined(__i386__) || defined(__x86_64__)
//#include <emmintrin.h>
//#endif
#include <inttypes.h>
#include <assert.h>

#include "macros.h"

#ifndef SPM_CRS_BITS
#define SPM_CRS_BITS 64
#endif

#include "macros.h"
#include "vector.h"
#include "spm_crs.h"
#include "mmf.h"
#include "spmv_method.h"

void *SPM_CRS_NAME(_init_mmf)(char *mmf_file,
                              uint64_t *nrows, uint64_t *ncols, uint64_t *nnz)
{

	SPM_CRS_TYPE *crs;
	crs = malloc(sizeof(SPM_CRS_TYPE));
	if (!crs){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	FILE *mmf = mmf_init(mmf_file, nrows, ncols, nnz);
	crs->nrows = *nrows;
	crs->ncols = *ncols;
	crs->nz = *nnz;

	/* allocate space for arrays */
	crs->values = malloc(sizeof(ELEM_TYPE)*crs->nz);
	crs->col_ind = malloc(sizeof(SPM_CRS_IDX_TYPE)*crs->nz);
	crs->row_ptr = malloc(sizeof(SPM_CRS_IDX_TYPE)*(crs->nrows + 1));
	if (!crs->values || !crs->col_ind || !crs->row_ptr){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	uint64_t row, col, row_prev;
	double val;
	uint64_t row_i=0, val_i=0;

	crs->row_ptr[row_i++] = val_i;
	row_prev = 0;
	while (mmf_get_next(mmf, &row, &col, &val)) {

		/* sanity checks */
		assert(row >= row_prev);
		assert(row < crs->nrows);
		assert(col < crs->ncols);
		assert(val_i < crs->nz);

		/* got a new row */
		if (row != row_prev){
			uint64_t i;
			/* handle empty rows, just in case */
			for (i=0; i<row - row_prev; i++){
				crs->row_ptr[row_i++] = val_i;
			}
			row_prev = row;
		}

		/* update values and colind arrays */
		crs->values[val_i] = (ELEM_TYPE)val;
		crs->col_ind[val_i] = (SPM_CRS_IDX_TYPE)col;
		val_i++;
	}
	crs->row_ptr[row_i++] = val_i;

	/* more sanity checks */
	assert(row_i == crs->nrows + 1);
	assert(val_i == crs->nz);
	fclose(mmf);

	return crs;
}

void SPM_CRS_NAME(_destroy)(void *spm)
{
	SPM_CRS_TYPE *crs = (SPM_CRS_TYPE *)spm;
	free(crs->values);
	free(crs->col_ind);
	free(crs->row_ptr);
	free(crs);
}

uint64_t SPM_CRS_NAME(_size)(void *spm)
{
	SPM_CRS_TYPE *crs = (SPM_CRS_TYPE *)spm;
	uint64_t ret = crs->nz*(sizeof(ELEM_TYPE) + sizeof(UINT_TYPE(SPM_CRS_BITS)));
	ret += (crs->nrows + 1)*sizeof(UINT_TYPE(SPM_CRS_BITS));
	return ret;
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

XSPMV_METH_INIT(
 SPM_CRS_NAME(_multiply),
 SPM_CRS_NAME(_init_mmf),
 SPM_CRS_NAME(_size),
 SPM_CRS_NAME(_destroy),
 sizeof(ELEM_TYPE)
)
