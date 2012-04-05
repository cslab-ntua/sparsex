/*
 * spm_vbl.c
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#ifndef SPM_CRS_BITS
#define SPM_CRS_BITS 64
#endif

#include "macros.h"
#include "vector.h"
#include "spm_crs.h"
#include "mmf.h"
#include "spmv_method.h"
#include "spm_vbl.h"

void *SPM_VBL_NAME(_init_mmf)(char *mmf_file,
                              uint64_t *nrows, uint64_t *ncols, uint64_t *nnz, void *metadata)
{

    SPM_CRS_TYPE *crs = SPM_CRS_NAME(_init_mmf)(mmf_file, nrows, ncols, nnz, metadata);
    SPM_VBL_TYPE *vbl = SPM_VBL_NAME(_init_crs)(crs);
    SPM_CRS_NAME(_destroy)(crs);
    return vbl;
}

void SPM_VBL_NAME(_destroy)(void *spm)
{
	SPM_VBL_TYPE *vbl = (SPM_VBL_TYPE *)spm;
	free(vbl->values);
	free(vbl->bcol_ind);
	free(vbl->row_ptr);
    free(vbl->bsize);
	free(vbl);
}

uint64_t SPM_VBL_NAME(_size)(void *spm)
{
	SPM_VBL_TYPE *vbl = (SPM_VBL_TYPE *) spm;
    uint64_t ret = vbl->nz * sizeof(*vbl->values) +
        vbl->nrows * sizeof(*vbl->row_ptr) +
        vbl->nblocks*(sizeof(*vbl->bcol_ind) + sizeof(*vbl->bsize));
	return ret;
}

void SPM_VBL_NAME(_multiply) (void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_VBL_TYPE *vbl = (SPM_VBL_TYPE *) spm;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *values = vbl->values;
	SPM_CRS_IDX_TYPE *row_ptr = vbl->row_ptr;
	SPM_CRS_IDX_TYPE *bcol_ind = vbl->bcol_ind;
    uint8_t *bsize = vbl->bsize;
	unsigned long n = vbl->nrows;
	register ELEM_TYPE yr;

	SPM_CRS_IDX_TYPE i, j, k, block_idx, blk_size;
    block_idx = 0;
    for (i = 0; i < n; i++) {
        yr = 0;
        for (j = row_ptr[i], blk_size = bsize[block_idx];
             j < row_ptr[i+1];
             j += blk_size, block_idx++, blk_size = bsize[block_idx]) {
            SPM_CRS_IDX_TYPE blk_start = bcol_ind[block_idx];
            for (k = 0; k < blk_size; k++)
                yr += values[j+k]*x[blk_start + k];
        }

        y[i] = yr;
    }

    return;
}

SPM_VBL_TYPE *SPM_VBL_NAME(_init_crs)(SPM_CRS_TYPE *crs)
{
	SPM_VBL_TYPE *vbl = malloc(sizeof(*vbl));
    if (!vbl) {
        fprintf(stderr, "malloc failed at %s:%d\n", __FILE__, __LINE__);
        exit(1);
    }
    vbl->nrows = crs->nrows;
    vbl->ncols = crs->ncols;
    vbl->nz = crs->nz;
    vbl->nblocks = 0;
    vbl->values = malloc(crs->nz*sizeof(*vbl->values));
    vbl->bcol_ind = malloc(crs->nz*sizeof(*vbl->bcol_ind));
    vbl->row_ptr = malloc(crs->nrows*sizeof(*vbl->row_ptr));
    vbl->bsize = malloc(crs->nz*sizeof(*vbl->bsize));

    SPM_CRS_IDX_TYPE i, j;
    SPM_CRS_IDX_TYPE col_prev, col_curr, col_delta;
    for (i = 0; i < crs->nrows; i++) {
        vbl->row_ptr[i]   = crs->row_ptr[i];
        vbl->row_ptr[i+1] = crs->row_ptr[i+1];

        /* That's a new row, start a new block. */
        col_curr = crs->col_ind[crs->row_ptr[i]];
        vbl->bcol_ind[vbl->nblocks] = col_curr;
        vbl->bsize[vbl->nblocks] = 1;

        for (j = crs->row_ptr[i], col_prev = crs->col_ind[j];
             j < crs->row_ptr[i+1]; j++) {
            vbl->values[j] = crs->values[j];
            col_curr = crs->col_ind[j];
            col_delta = col_curr - col_prev;
            if (col_delta == 0)
                continue;
            else if (col_delta == 1) {
                int new_block_size = vbl->bsize[vbl->nblocks] + 1;
                if (new_block_size > UINT8_MAX) {
                    /* Finalize the current block ... */
                    vbl->bsize[vbl->nblocks] = UINT8_MAX;
                    /* ... and start a new one */
                    vbl->nblocks++;
                    vbl->bcol_ind[vbl->nblocks] = col_curr;
                    vbl->bsize[vbl->nblocks] = (uint8_t)
                        (new_block_size - UINT8_MAX);
                } else
                    vbl->bsize[vbl->nblocks] = new_block_size;
            } else {
                /* That's a new block */
                vbl->nblocks++;
                vbl->bcol_ind[vbl->nblocks] = col_curr;
                vbl->bsize[vbl->nblocks] = 1;
            }

            col_prev = col_curr;
        }

        /* Line finished, increment block counter */
        vbl->nblocks++;
        assert(vbl->nblocks <= crs->nz);
    }

    /* Truncate the extra space for the bcol_ind and bsize array. */
    vbl->bcol_ind = realloc(vbl->bcol_ind,
                            vbl->nblocks*sizeof(*vbl->bcol_ind));
    vbl->bsize = realloc(vbl->bsize, (vbl->nblocks+1)*sizeof(*vbl->bsize));
    vbl->bsize[vbl->nblocks] = 0;
    return vbl;
}


XSPMV_METH_INIT(
 SPM_VBL_NAME(_multiply),
 SPM_VBL_NAME(_init_mmf),
 SPM_VBL_NAME(_size),
 SPM_VBL_NAME(_destroy),
 sizeof(ELEM_TYPE)
)
