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

#ifdef SPM_NUMA

#include <numa.h>
#include "numa_util.h"

void *SPM_BCSR_MT_NAME(_numa_init_mmf)(char *mmf_file,
                                      uint64_t *rows_nr, uint64_t *cols_nr,
                                      uint64_t *nz_nr, void *metadata)
{

    spm_mt_t *spm_mt = SPM_BCSR_MT_NAME(_init_mmf)(mmf_file,
                                                  rows_nr, cols_nr,
                                                  nz_nr, metadata);
    int nr_threads = spm_mt->nr_threads;
    size_t *bvalues_parts = malloc(nr_threads*sizeof(*bvalues_parts));
    size_t *browptr_parts = malloc(nr_threads*sizeof(*browptr_parts));
    size_t *bcolind_parts = malloc(nr_threads*sizeof(*bcolind_parts));
    int *nodes = malloc(nr_threads*sizeof(*nodes));

    // just reallocate in a numa-aware fashion the data structures
    int i;
    SPM_BCSR_TYPE *bcsr = NULL;
    for (i = 0; i < nr_threads; i++) {
        spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
        SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *) spm_thread->spm;
        bcsr = bcsr_mt->bcsr;
        SPM_CRS_IDX_TYPE row_start = bcsr_mt->row_start;
        SPM_CRS_IDX_TYPE row_end = bcsr_mt->row_end;
        SPM_CRS_IDX_TYPE brow_start = row_start / bcsr->br;
        SPM_CRS_IDX_TYPE brow_end = row_end / bcsr->br;
        SPM_CRS_IDX_TYPE vstart = bcsr->brow_ptr[brow_start];
        SPM_CRS_IDX_TYPE vend = bcsr->brow_ptr[brow_end];
        SPM_CRS_IDX_TYPE nblocks_part = (vend - vstart) / (bcsr->br*bcsr->bc);
        browptr_parts[i] = (brow_end-brow_start)*sizeof(*bcsr->brow_ptr);
        bcolind_parts[i] = nblocks_part*sizeof(*bcsr->bcol_ind);
        bvalues_parts[i] = (vend-vstart)*sizeof(*bcsr->bvalues);
        nodes[i] = numa_node_of_cpu(spm_thread->cpu);
        spm_thread->node = nodes[i];
        spm_thread->row_start = row_start;
        spm_thread->nr_rows = row_end - row_start;
    }

    // sanity check (+ get rid of compiler warning about uninit. variable)
    assert(bcsr);

    SPM_CRS_IDX_TYPE nr_values = bcsr->nblocks*bcsr->br*bcsr->bc;
    SPM_CRS_IDX_TYPE *new_browptr =
        alloc_interleaved((bcsr->nbrows+1)*sizeof(*bcsr->brow_ptr),
                          browptr_parts, nr_threads, nodes);

    SPM_CRS_IDX_TYPE *new_bcolind =
        alloc_interleaved(bcsr->nblocks*sizeof(*bcsr->bcol_ind),
                          bcolind_parts, nr_threads, nodes);
    ELEM_TYPE *new_bvalues =
        alloc_interleaved(nr_values*sizeof(*bcsr->bvalues),
                          bvalues_parts, nr_threads, nodes);

    // copy old data to the new one
    memcpy(new_browptr, bcsr->brow_ptr,
           (bcsr->nbrows+1)*sizeof(*bcsr->brow_ptr));
    memcpy(new_bcolind, bcsr->bcol_ind, bcsr->nblocks*sizeof(*bcsr->bcol_ind));
    memcpy(new_bvalues, bcsr->bvalues, nr_values*sizeof(*bcsr->bvalues));

    // check allocation
	int alloc_err;
	alloc_err = check_interleaved((void *) new_browptr, browptr_parts,
	                              nr_threads, nodes);
	print_alloc_status("BCSR browptr", alloc_err);
	alloc_err = check_interleaved((void *) new_bcolind,  bcolind_parts,
	                              nr_threads, nodes);
	print_alloc_status("BCSR bcolind", alloc_err);
	alloc_err = check_interleaved((void *) new_bvalues, bvalues_parts,
	                              nr_threads, nodes);
	print_alloc_status("BCSR bvalues", alloc_err);

    // free old data and replace with the new one
    free(bcsr->brow_ptr);
    free(bcsr->bcol_ind);
    free(bcsr->bvalues);
    bcsr->brow_ptr = new_browptr;
    bcsr->bcol_ind = new_bcolind;
    bcsr->bvalues = new_bvalues;

    // free the auxiliaries
    free(browptr_parts);
    free(bcolind_parts);
    free(bvalues_parts);
    free(nodes);
	return spm_mt;
}

void SPM_BCSR_MT_NAME(_numa_destroy)(void *spm)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_BCSR_MT_TYPE *bcsr_mt = (SPM_BCSR_MT_TYPE *) spm_thread->spm;
    SPM_BCSR_TYPE *bcsr = bcsr_mt->bcsr;
    SPM_CRS_IDX_TYPE nr_values = bcsr->nblocks*bcsr->br*bcsr->bc;
    free_interleaved(bcsr->brow_ptr, (bcsr->nbrows+1)*sizeof(*bcsr->brow_ptr));
    free_interleaved(bcsr->bcol_ind, bcsr->nblocks*sizeof(*bcsr->bcol_ind));
    free_interleaved(bcsr->bvalues, nr_values*sizeof(*bcsr->bvalues));
    free(bcsr);
    free(bcsr_mt);
    free(spm_thread);
    free(spm_mt);
}

uint64_t SPM_BCSR_MT_NAME(_numa_size)(void *spm)
{
    return SPM_BCSR_MT_NAME(_size)(spm);
}

void SPM_BCSR_MT_NAME(_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
    SPM_BCSR_MT_NAME(_multiply)(spm, in, out);
}

XSPMV_MT_METH_INIT(
 SPM_BCSR_MT_NAME(_numa_multiply),
 SPM_BCSR_MT_NAME(_numa_init_mmf),
 SPM_BCSR_MT_NAME(_numa_size),
 SPM_BCSR_MT_NAME(_numa_destroy),
 sizeof(ELEM_TYPE)
)

#endif /* SPM_NUMA */
