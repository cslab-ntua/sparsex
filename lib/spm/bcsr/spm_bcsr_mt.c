/*
 * spm_bcsr_mt.c -- Multithreaded BCSR
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <assert.h>
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
#include "spm_bcsr.h"
#include "spm_bcsr_mt.h"
#include "mult/multiply.h"

static void SPM_BCSR_MT_NAME(_split)(SPM_BCSR_TYPE *mat,
                                     SPM_CRS_IDX_TYPE nr_splits,
                                     SPM_CRS_IDX_TYPE *splits);

void *SPM_BCSR_MT_NAME(_init_mmf)(char *mmf_file,
                                 uint64_t *rows_nr, uint64_t *cols_nr,
                                 uint64_t *nz_nr, void *metadata)
{
	int i;
	unsigned int nr_cpus, *cpus;
	spm_mt_t *spm_mt;
	spm_mt_thread_t *spm_thread;
	SPM_BCSR_TYPE *bcsr;
	SPM_BCSR_MT_TYPE *bcsr_mt;
	SPM_CRS_TYPE *crs;

    assert(metadata);

    // set affinity of the current thread
	mt_get_options(&nr_cpus, &cpus);
    setaffinity_oncpu(cpus[0]);

	crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr, metadata);
    bcsr_metadata_t *meta = (bcsr_metadata_t *) metadata;
    bcsr = SPM_BCSR_NAME(_init_crs)(crs, meta->br, meta->bc,
                                    BCSR_PAD_MAX(meta->br, meta->bc));
    
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

	bcsr_mt = malloc(sizeof(SPM_BCSR_MT_TYPE)*nr_cpus);
	if ( !bcsr_mt ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

    SPM_CRS_IDX_TYPE *splits = malloc((nr_cpus+1)*sizeof(*splits));

    SPM_BCSR_MT_NAME(_split)(bcsr, nr_cpus, splits);
    for (i = 0; i < nr_cpus; i++) {
        bcsr_mt[i].row_start = splits[i];
        bcsr_mt[i].row_end = splits[i+1];
/*         bcsr_mt[i].nnz_nr = */
/*             bcsr->brow_ptr[splits[i+1]] - bcsr->brow_ptr[splits[i]]; */
        bcsr_mt[i].bcsr = bcsr;
        spm_thread = spm_mt->spm_threads + i;
        spm_thread->cpu = cpus[i];
        spm_thread->spm = bcsr_mt + i;
    }

    free(splits);
	free(cpus);
    SPM_CRS_NAME(_destroy)(crs);
	return spm_mt;
}

void SPM_BCSR_MT_NAME(_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *) spm_thread->spm;
	SPM_BCSR_NAME(_destroy)(bcsr_mt->bcsr);
	free(bcsr_mt);
	free(spm_thread);
	free(spm_mt);
}

uint64_t SPM_BCSR_MT_NAME(_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *)spm_thread->spm;
	return SPM_BCSR_NAME(_size)(bcsr_mt->bcsr);
}

void SPM_BCSR_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *) spm;
    SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
    SPM_BCSR_MT_NAME(_mult_method)(bcsr->br, bcsr->bc)(spm, in, out);    
}

XSPMV_MT_METH_INIT(
 SPM_BCSR_MT_NAME(_multiply),
 SPM_BCSR_MT_NAME(_init_mmf),
 SPM_BCSR_MT_NAME(_size),
 SPM_BCSR_MT_NAME(_destroy),
 sizeof(ELEM_TYPE)
)

static void SPM_BCSR_MT_NAME(_split)(SPM_BCSR_TYPE *mat,
                                     SPM_CRS_IDX_TYPE nr_splits,
                                     SPM_CRS_IDX_TYPE *splits)
{
    SPM_CRS_IDX_TYPE nr_nzeros = mat->nblocks*mat->br*mat->bc;
    SPM_CRS_IDX_TYPE *brow_ptr = mat->brow_ptr;
    SPM_CRS_IDX_TYPE min_nzeros_per_split = nr_nzeros / nr_splits;
    SPM_CRS_IDX_TYPE i;
    SPM_CRS_IDX_TYPE nr_nzeros_per_split = 0;
    SPM_CRS_IDX_TYPE split_pos = 1;

    splits[0] = 0;
    /* Iterate over block rows */
    for (i = 0; i < mat->nbrows; i++) {
        nr_nzeros_per_split += brow_ptr[i+1] - brow_ptr[i];
        if (nr_nzeros_per_split > min_nzeros_per_split) {
            /* new split point found */
            splits[split_pos++] = (i+1)*mat->br;
            nr_nzeros_per_split -= min_nzeros_per_split;
        }
    }

    splits[nr_splits] = mat->nbrows * mat->br;
    return;
}
