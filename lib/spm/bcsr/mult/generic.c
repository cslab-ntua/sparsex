/*
 * bcsr/generic.c -- BCSR generic multiplication method
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "mult/multiply.h"

void SPM_BCSR_NAME(_multiply_generic) (void *spm, VECTOR_TYPE *in,
                                       VECTOR_TYPE *out)
{
	SPM_BCSR_TYPE *mat = (SPM_BCSR_TYPE *) spm;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *bvalues = mat->bvalues;
	SPM_CRS_IDX_TYPE *brow_ptr = mat->brow_ptr;
	SPM_CRS_IDX_TYPE *bcol_ind = mat->bcol_ind;
	uint64_t r = mat->br;
	uint64_t c = mat->bc;
	const SPM_CRS_IDX_TYPE r_start = 0;
	const SPM_CRS_IDX_TYPE r_end = mat->nrows;

	SPM_CRS_IDX_TYPE i, _i, j, _j, k, l;
	register ELEM_TYPE yr;
	for (i = r_start, _i = r_start / r; i < r_end; i += r, _i++) {
		for (j = brow_ptr[_i], _j = j / (r*c); j < brow_ptr[_i+1];
			 j += r*c, _j++) {
			SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
			for (k = 0; k < r; k++) {
				yr = 0;
				for (l = 0; l < c; l++)
					yr += bvalues[j+k*c+l]*x[x_start + l];
				y[i+k] += yr;
			}
		}
	}

	return;
}

void SPM_BCSR_MT_NAME(_multiply_generic) (void *spm, VECTOR_TYPE *in,
                                          VECTOR_TYPE *out)
{
	SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *) spm;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *bvalues = bcsr_mt->bcsr->bvalues;
	SPM_CRS_IDX_TYPE *brow_ptr = bcsr_mt->bcsr->brow_ptr;
	SPM_CRS_IDX_TYPE *bcol_ind = bcsr_mt->bcsr->bcol_ind;
	uint64_t r = bcsr_mt->bcsr->br;
	uint64_t c = bcsr_mt->bcsr->bc;
	const SPM_CRS_IDX_TYPE r_start = bcsr_mt->row_start;
	const SPM_CRS_IDX_TYPE r_end = bcsr_mt->row_end;

	SPM_CRS_IDX_TYPE i, _i, j, _j, k, l;
	register ELEM_TYPE yr;
	for (i = r_start, _i = r_start / r; i < r_end; i += r, _i++) {
		for (j = brow_ptr[_i], _j = j / (r*c); j < brow_ptr[_i+1];
			 j += r*c, _j++) {
			SPM_CRS_IDX_TYPE x_start = bcol_ind[_j];
			for (k = 0; k < r; k++) {
				yr = 0;
				for (l = 0; l < c; l++)
					yr += bvalues[j+k*c+l]*x[x_start + l];
				y[i+k] += yr;
			}
		}
	}

	return;
}
