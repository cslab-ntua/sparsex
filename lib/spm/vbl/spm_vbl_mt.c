/*
 * spm_vbl_mt.c -- Multithreaded VBL
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
    SPM_CRS_IDX_TYPE *bsplits = malloc((nr_cpus+1)*sizeof(*bsplits));

    SPM_VBL_MT_NAME(_split)(vbl, nr_cpus, splits, bsplits);
    for (i = 0; i < nr_cpus; i++) {
        vbl_mt[i].row_start = splits[i];
        vbl_mt[i].row_end = splits[i+1];
        vbl_mt[i].bstart = bsplits[i];
        vbl_mt[i].bend = bsplits[i+1];
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
    blk_splits[nr_splits] = vbl->nblocks;

//	for (i = 0; i <= nr_splits; i++)
//		printf("splits[%d] = %ld\n", i, splits[i]);
    return;
}

#ifdef SPM_NUMA
#include <numa.h>
#include "numa_util.h"

void *SPM_VBL_MT_NAME(_numa_init_mmf)(char *mmf_file,
                                      uint64_t *rows_nr, uint64_t *cols_nr,
                                      uint64_t *nz_nr, void *metadata)
{
    spm_mt_t *spm_mt = SPM_VBL_MT_NAME(_init_mmf)(mmf_file,
                                                  rows_nr, cols_nr,
                                                  nz_nr, metadata);

    int nr_threads = spm_mt->nr_threads;
    size_t *values_parts = malloc(nr_threads*sizeof(*values_parts));
    size_t *rowptr_parts = malloc(nr_threads*sizeof(*rowptr_parts));
    size_t *bcolind_parts = malloc(nr_threads*sizeof(*bcolind_parts));
    size_t *bsize_parts = malloc(nr_threads*sizeof(*bsize_parts));
    int *nodes = malloc(nr_threads*sizeof(*nodes));

    // just reallocate in a numa-aware fashion the data structures
    int i;
    SPM_VBL_TYPE *vbl = NULL;
    for (i = 0; i < nr_threads; i++) {
        spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
        SPM_VBL_MT_TYPE *vbl_mt = (SPM_VBL_MT_TYPE *) spm_thread->spm;
        vbl = vbl_mt->vbl;
        SPM_CRS_IDX_TYPE row_start = vbl_mt->row_start;
        SPM_CRS_IDX_TYPE row_end = vbl_mt->row_end;
        SPM_CRS_IDX_TYPE bstart = vbl_mt->bstart;
        SPM_CRS_IDX_TYPE bend = vbl_mt->bend;
        SPM_CRS_IDX_TYPE vstart = vbl->row_ptr[vbl_mt->row_start];
        SPM_CRS_IDX_TYPE vend = vbl->row_ptr[vbl_mt->row_end];
        rowptr_parts[i] = (row_end-row_start)*sizeof(*vbl->row_ptr);
        bcolind_parts[i] = (bend-bstart)*sizeof(*vbl->bcol_ind);
        bsize_parts[i] = (bend-bstart)*sizeof(*vbl->bsize);
        values_parts[i] = (vend-vstart)*sizeof(*vbl->values);
        nodes[i] = numa_node_of_cpu(spm_thread->cpu);
        spm_thread->node = nodes[i];
        spm_thread->row_start = row_start;
        spm_thread->nr_rows = row_end - row_start;
    }

    // sanity check (+ get rid of compiler warning about uninit. variable)
    assert(vbl);

    SPM_CRS_IDX_TYPE *new_rowptr =
        alloc_interleaved((vbl->nrows+1)*sizeof(*vbl->row_ptr),
                          rowptr_parts, nr_threads, nodes);

    SPM_CRS_IDX_TYPE *new_bcolind =
        alloc_interleaved(vbl->nblocks*sizeof(*vbl->bcol_ind),
                          bcolind_parts, nr_threads, nodes);
    uint8_t *new_bsize = alloc_interleaved(vbl->nblocks*sizeof(*vbl->bsize),
                                           bsize_parts, nr_threads, nodes);
    ELEM_TYPE *new_values = alloc_interleaved(vbl->nz*sizeof(*vbl->values),
                                              values_parts, nr_threads,
                                              nodes);

    // copy old data to the new one
    memcpy(new_rowptr, vbl->row_ptr, (vbl->nrows+1)*sizeof(*vbl->row_ptr));
    memcpy(new_bcolind, vbl->bcol_ind, vbl->nblocks*sizeof(*vbl->bcol_ind));
    memcpy(new_bsize, vbl->bsize, vbl->nblocks*sizeof(*vbl->bsize));
    memcpy(new_values, vbl->values, vbl->nz*sizeof(*vbl->values));

    // check allocation
	int alloc_err;
	alloc_err = check_interleaved((void *) new_rowptr, rowptr_parts,
	                              nr_threads, nodes);
	print_alloc_status("VBL rowptr", alloc_err);
	alloc_err = check_interleaved((void *) new_bcolind,  bcolind_parts,
	                              nr_threads, nodes);
	print_alloc_status("VBL bcolind", alloc_err);
	alloc_err = check_interleaved((void *) new_values, values_parts,
	                              nr_threads, nodes);
	print_alloc_status("VBL values", alloc_err);

    // free old data and replace with the new one
    free(vbl->row_ptr);
    free(vbl->bcol_ind);
    free(vbl->bsize);
    free(vbl->values);
    vbl->row_ptr = new_rowptr;
    vbl->bcol_ind = new_bcolind;
    vbl->bsize = new_bsize;
    vbl->values = new_values;

    // free the auxiliaries
    free(rowptr_parts);
    free(bcolind_parts);
    free(bsize_parts);
    free(values_parts);
    free(nodes);
    return spm_mt;
}

void SPM_VBL_MT_NAME(_numa_destroy)(void *spm)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_VBL_MT_TYPE *vbl_mt = (SPM_VBL_MT_TYPE *) spm_thread->spm;
    SPM_VBL_TYPE *vbl = vbl_mt->vbl;
    free_interleaved(vbl->row_ptr, (vbl->nrows+1)*sizeof(*vbl->row_ptr));
    free_interleaved(vbl->bcol_ind, vbl->nblocks*sizeof(*vbl->bcol_ind));
    free_interleaved(vbl->bsize, vbl->nblocks*sizeof(*vbl->bsize));
    free_interleaved(vbl->values, vbl->nz*sizeof(*vbl->values));
    free(vbl);
    free(vbl_mt);
    free(spm_thread);
    free(spm_mt);
}

uint64_t SPM_VBL_MT_NAME(_numa_size)(void *spm)
{
    return SPM_VBL_MT_NAME(_size)(spm);
}

void SPM_VBL_MT_NAME(_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
    SPM_VBL_MT_NAME(_multiply)(spm, in, out);
}

XSPMV_MT_METH_INIT(
 SPM_VBL_MT_NAME(_numa_multiply),
 SPM_VBL_MT_NAME(_numa_init_mmf),
 SPM_VBL_MT_NAME(_numa_size),
 SPM_VBL_MT_NAME(_numa_destroy),
 sizeof(ELEM_TYPE)
)

#endif  /* SPM_NUMA */
