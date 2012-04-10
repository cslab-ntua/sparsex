/*
 * spm_bcsr.c -- BCSR implemetation
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
#include "spm_bcsr.h"
#include "blocks.h"
#include "bitstr.h"
#include "util.h"
#include "mult/multiply.h"

static void _scan_csr(const SPM_CRS_TYPE *mat,
                      SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c,
                      SPM_CRS_IDX_TYPE max_pad, SPM_CRS_IDX_TYPE *nr_blocks,
                      SPM_CRS_IDX_TYPE *nr_leftovers, block_cache_t *blk_cache,
                      int skip);

static void _fill_bcsr(SPM_BCSR_TYPE *mat, const SPM_CRS_TYPE *matin,
                       SPM_CRS_IDX_TYPE max_pad,
                       SPM_CRS_TYPE *leftover, block_cache_t *blk_cache);

void *SPM_BCSR_NAME(_init_mmf)(char *mmf_file,
                               uint64_t *nrows, uint64_t *ncols, uint64_t *nnz,
                               void *metadata)
{

	assert(metadata);
	SPM_CRS_TYPE *crs = SPM_CRS_NAME(_init_mmf)(mmf_file, nrows, ncols, nnz,
	                                            NULL);
	bcsr_metadata_t *meta = (bcsr_metadata_t *) metadata;
	SPM_BCSR_TYPE *bcsr =
	    SPM_BCSR_NAME(_init_crs)(crs, meta->br, meta->bc,
	                             BCSR_PAD_MAX(meta->br, meta->bc));
	SPM_CRS_NAME(_destroy)(crs);
	return bcsr;
}

void SPM_BCSR_NAME(_destroy)(void *spm)
{
	SPM_BCSR_TYPE *bcsr = (SPM_BCSR_TYPE *) spm;
	free(bcsr->bvalues);
	free(bcsr->bcol_ind);
	free(bcsr->brow_ptr);
	free(bcsr);
}

uint64_t SPM_BCSR_NAME(_size)(void *spm)
{
	SPM_BCSR_TYPE *mat = (SPM_BCSR_TYPE *) spm;
//	  SPM_CRS_IDX_TYPE nbrows = iceil(mat->nrows, mat->br);
	uint64_t size =
	    mat->br * mat->bc * mat->nblocks * sizeof(*mat->bvalues) +
	    mat->nbrows * sizeof(*mat->brow_ptr) +
	    mat->nblocks * sizeof(*mat->bcol_ind);
	return size;
}

void SPM_BCSR_NAME(_multiply) (void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_BCSR_TYPE *mat = (SPM_BCSR_TYPE *) spm;
	SPM_BCSR_NAME(_mult_method)(mat->br, mat->bc)(spm, in, out);
}

SPM_BCSR_TYPE *SPM_BCSR_NAME(_init_crs)(const SPM_CRS_TYPE *crs,
                                        uint64_t r, uint64_t c, uint64_t pad)
{
	SPM_CRS_IDX_TYPE nr_blocks, nr_leftovers;
	block_cache_t	*blk_cache;
	SPM_CRS_IDX_TYPE blk_cache_size;
	blk_cache_size = iceil(crs->ncols, c);
	blk_cache      = block_cache_create(blk_cache_size, r, c,
	                                    BLOCK_STORAGE_RW);

	/* Scan the input matrix for blocks */
	_scan_csr(crs, r, c, pad, &nr_blocks, &nr_leftovers, blk_cache,
	          pad == BCSR_PAD_MAX(r,c));

	assert(!nr_leftovers && pad == BCSR_PAD_MAX(r,c));
	SPM_BCSR_TYPE *bcsr = malloc(sizeof(*bcsr));
	bcsr->nbrows = iceil(crs->nrows, r);
	bcsr->brow_ptr = malloc((bcsr->nbrows+1)*sizeof(*bcsr->brow_ptr));
	bcsr->bcol_ind = malloc(nr_blocks*sizeof(*bcsr->bcol_ind));
	bcsr->bvalues = malloc(nr_blocks*r*c*sizeof(*bcsr->bvalues));
	bcsr->nz = crs->nz;
	bcsr->nblocks = nr_blocks;
	bcsr->nrows = crs->nrows;
	bcsr->ncols = crs->ncols;
	bcsr->br = r;
	bcsr->bc = c;
	bcsr->storage = blk_layout_list[blk_layout_default].sl_layout;

	_fill_bcsr(bcsr, crs, pad, NULL, blk_cache);

	block_cache_delete(blk_cache);
	return bcsr;
}

XSPMV_METH_INIT(
	SPM_BCSR_NAME(_multiply),
	SPM_BCSR_NAME(_init_mmf),
	SPM_BCSR_NAME(_size),
	SPM_BCSR_NAME(_destroy),
	sizeof(ELEM_TYPE)
)

static void
_scan_csr(const SPM_CRS_TYPE *crs, SPM_CRS_IDX_TYPE r, SPM_CRS_IDX_TYPE c,
          SPM_CRS_IDX_TYPE max_pad, SPM_CRS_IDX_TYPE *nr_blocks,
          SPM_CRS_IDX_TYPE *nr_leftovers, block_cache_t *blk_cache, int skip)
{
	SPM_CRS_IDX_TYPE i, j, k;
	spm_elem_t curr_elem;

	*nr_blocks	  = 0;
	*nr_leftovers = 0;
	for (i = 0; i < crs->nrows; i += r) {
		for (j = 0; j < r; j++) {
			curr_elem.row = i + j;
			if (curr_elem.row + 1 > crs->nrows)
				break;
			for (k = crs->row_ptr[i + j];
				 k < crs->row_ptr[i + j + 1]; k++) {
				curr_elem.col	= crs->col_ind[k];
				curr_elem.value = crs->values[k];
				block_cache_add_elem(blk_cache, &curr_elem);
			}
		}

		*nr_blocks += block_cache_nr_blocks(blk_cache);

		block_cache_iter_t	iter;
		block_t				block;

		/* Skip the search for leftovers */
		if (!skip) {
			block_cache_iter_init(&iter, blk_cache);
			block_cache_iter_begin(&iter);
			while (block_cache_iter_more(&iter)) {
				block_cache_iter_next(&iter, &block);
				if (block.pad > max_pad) {
					*nr_leftovers += r*c - block.pad;
					(*nr_blocks)--;
				}
			}
		}

		block_cache_clear(blk_cache);	/* clear the cache */
	}

	return;
}

/*
 *	FIXME: The filling algorithm is very similar to the scan algorithm. Maybe
 *	merge to avoid code replication or use the template pattern?
 */
static void
_fill_bcsr(SPM_BCSR_TYPE *mat, const SPM_CRS_TYPE *matin,
           SPM_CRS_IDX_TYPE max_pad,
           SPM_CRS_TYPE *leftover, block_cache_t *blk_cache)
{
	SPM_CRS_IDX_TYPE	i, j, k;
	spm_elem_t curr_elem;
	SPM_CRS_IDX_TYPE	r;
	bcsr_context_t	bcsr_context;
	csr_context_t	csr_context;

	bcsr_context.next_brow	= 0;
	bcsr_context.next_block = 0;
	csr_context.next_row	= 0;
	csr_context.next_val	= 0;
	csr_context.cache		= blk_cache;

	r = mat->br;
	for (i = 0; i < matin->nrows; i += r) {
		for (j = 0; j < r; j++) {
			curr_elem.row = i + j;
			if (curr_elem.row + 1 > matin->nrows)
				break;
			for (k = matin->row_ptr[i + j];
				 k < matin->row_ptr[i + j + 1]; k++) {
				curr_elem.col	= matin->col_ind[k];
				curr_elem.value = matin->values[k];
				block_cache_add_elem(blk_cache, &curr_elem);
			}
		}

		block_cache_iter_t	iter;
		block_t             block;

		/*
		 * First sort the valid blocks. This is not a necessary step; just to
		 * be strictly conformant to the BCSR format.
		 */
		block_cache_sort(blk_cache);
		block_cache_iter_init(&iter, blk_cache);
		block_cache_iter_begin(&iter);
		while (block_cache_iter_more(&iter)) {
			block_cache_iter_next(&iter, &block);
			if (block.pad > max_pad) {
				/* This block is added to the leftover */
				if (leftover)
					csr_add(leftover, &block, &csr_context);
			} else {
				bcsr_add(mat, &block, &bcsr_context);
			}
		}

		block_cache_clear(blk_cache);	/* clear the cache */
	}

	bcsr_finalize(mat, &bcsr_context);
	if (leftover)
		csr_finalize(leftover, &csr_context);

	return;
}
