/*
 * spm_vbl_mt.c -- Multithreaded VBL
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdint.h>

#include "macros.h"
#include "mt_lib.h"
#include "spm_mt.h"
#include "spm_crs.h"
#include "spm_crs_mt.h"
#include "spmv_method.h"
#include "spm_vbl.h"
#include "spm_vbl_mt.h"

static void SPM_VBL_MT_NAME(_split)(SPM_VBL_TYPE *vbl,
                                    SPM_CRS_IDX_TYPE nr_splits,
                                    SPM_CRS_IDX_TYPE *splits,
                                    SPM_CRS_IDX_TYPE *blk_splits);

void *SPM_VBL_MT_NAME(_init_mmf)(char *mmf_file,
                                 uint64_t *rows_nr, uint64_t *cols_nr,
                                 uint64_t *nz_nr, void *metadata)
{
	int i;
	unsigned int nr_cpus, *cpus;
	spm_mt_t *spm_mt;
	spm_mt_thread_t *spm_thread;
	SPM_VBL_TYPE *vbl;
	SPM_VBL_MT_TYPE *vbl_mt;
	SPM_CRS_TYPE *crs;

    // set affinity of the current thread
	mt_get_options(&nr_cpus, &cpus);
    setaffinity_oncpu(cpus[0]);

	crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr, metadata);
    vbl = SPM_VBL_NAME(_init_crs)(crs);
    
	spm_mt = malloc(sizeof(spm_mt_t));
	if ( !spm_mt ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	spm_mt->nr_threads = nr_cpus;
	spm_mt->spm_threads = malloc(sizeof(spm_mt_thread_t)*nr_cpus);
	if ( !spm_mt->spm_threads ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	vbl_mt = malloc(sizeof(SPM_VBL_MT_TYPE)*nr_cpus);
	if ( !vbl_mt ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

    SPM_CRS_IDX_TYPE *splits = malloc((nr_cpus+1)*sizeof(*splits));
    SPM_CRS_IDX_TYPE *bsplits = malloc(nr_cpus*sizeof(*bsplits));

    SPM_VBL_MT_NAME(_split)(vbl, nr_cpus, splits, bsplits);
    for (i = 0; i < nr_cpus; i++) {
        vbl_mt[i].row_start = splits[i];
        vbl_mt[i].row_end = splits[i+1];
        vbl_mt[i].bstart = bsplits[i];
        vbl_mt[i].nnz_nr =
            vbl->row_ptr[splits[i+1]] - vbl->row_ptr[splits[i]];
        vbl_mt[i].vbl = vbl;
        spm_thread = spm_mt->spm_threads + i;
        spm_thread->cpu = cpus[i];
        spm_thread->spm = vbl_mt + i;
    }

    free(splits);
    free(bsplits);
	free(cpus);
    SPM_CRS_NAME(_destroy)(crs);
	return spm_mt;
}

void SPM_VBL_MT_NAME(_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_VBL_MT_TYPE *vbl_mt = (SPM_VBL_MT_TYPE *)spm_thread->spm;
	SPM_VBL_NAME(_destroy)(vbl_mt->vbl);
	free(vbl_mt);
	free(spm_thread);
	free(spm_mt);
}

uint64_t SPM_VBL_MT_NAME(_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_VBL_MT_TYPE *vbl_mt = (SPM_VBL_MT_TYPE *)spm_thread->spm;
	return SPM_VBL_NAME(_size)(vbl_mt->vbl);
}

void SPM_VBL_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_VBL_MT_TYPE *vbl_mt = (SPM_VBL_MT_TYPE *) spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = vbl_mt->vbl->values;
	SPM_CRS_IDX_TYPE *row_ptr = vbl_mt->vbl->row_ptr;
	SPM_CRS_IDX_TYPE *bcol_ind = vbl_mt->vbl->bcol_ind;
    uint8_t *bsize = vbl_mt->vbl->bsize;
	const unsigned long row_start = vbl_mt->row_start;
	const unsigned long row_end = vbl_mt->row_end;
    const unsigned long bstart = vbl_mt->bstart;
	register ELEM_TYPE yr;

	unsigned long i, j, k, block_idx, blk_size;

    block_idx = bstart;
    for (i = row_start; i < row_end; i++) {
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
}
XSPMV_MT_METH_INIT(
 SPM_VBL_MT_NAME(_multiply),
 SPM_VBL_MT_NAME(_init_mmf),
 SPM_VBL_MT_NAME(_size),
 SPM_VBL_MT_NAME(_destroy),
 sizeof(ELEM_TYPE)
)

static void SPM_VBL_MT_NAME(_split)(SPM_VBL_TYPE *vbl,
                                    SPM_CRS_IDX_TYPE nr_splits,
                                    SPM_CRS_IDX_TYPE *splits,
                                    SPM_CRS_IDX_TYPE *blk_splits)
{
    SPM_CRS_IDX_TYPE nr_nzeros = vbl->nz;
    SPM_CRS_IDX_TYPE *row_ptr = vbl->row_ptr;
    SPM_CRS_IDX_TYPE min_nzeros_per_split = nr_nzeros / nr_splits;
    SPM_CRS_IDX_TYPE i, j;

    SPM_CRS_IDX_TYPE nr_nzeros_per_split = 0;
    SPM_CRS_IDX_TYPE split_pos = 1;
    SPM_CRS_IDX_TYPE blk_idx = 0;
    SPM_CRS_IDX_TYPE blk_size;

    splits[0]     = 0;
    blk_splits[0] = 0;
    for (i = 0; i < vbl->nrows; i++) {
        nr_nzeros_per_split += row_ptr[i+1] - row_ptr[i];

        /* Find the number of blocks in this row */
        for (j = row_ptr[i], blk_size = vbl->bsize[blk_idx];
             j < row_ptr[i+1];
             j += blk_size, blk_idx++, blk_size = vbl->bsize[blk_idx])
            ;

        if (nr_nzeros_per_split > min_nzeros_per_split) {
            /* new split point found */
            splits[split_pos] = i + 1;
            blk_splits[split_pos] = blk_idx;
            split_pos++;
            nr_nzeros_per_split -= min_nzeros_per_split;
        }
    }

    splits[nr_splits] = vbl->nrows;

//	for (i = 0; i <= nr_splits; i++)
//		printf("splits[%d] = %ld\n", i, splits[i]);
    return;
}
