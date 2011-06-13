/*
 * spm_crs_sym.c -- Implementation of expansion of CSR for symmetric sparse matrices.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "spm_crs_sym_mt.h"

void *SPM_CRS_SYM_MT_NAME(_init_mmf)(char *mmf_file,
                                     uint64_t *nrows, uint64_t *ncols,
                                     uint64_t *nnz)
{
    unsigned int i;
    uint32_t ncpus, *cpus;
    uint64_t cur_row, row_start, elems_limit, elems_total;
    spm_mt_t *spm_mt;
    spm_mt_thread_t *spm_thread;
    SPM_CRS_SYM_MT_TYPE *crs_mt;
    SPM_CRS_SYM_TYPE *crs;
    
    crs = (SPM_CRS_SYM_TYPE *) SPM_CRS_SYM_NAME(_init_mmf)(mmf_file, nrows, ncols, nnz);
    mt_get_options(&ncpus, &cpus);
    
    spm_mt = (spm_mt_t *) malloc(sizeof(spm_mt_t));
    if (!spm_mt) {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }
    
    spm_mt->nr_threads = ncpus;
    spm_mt->spm_threads = (spm_mt_thread_t *) malloc(ncpus * sizeof(spm_mt_thread_t));
    if (!spm_mt->spm_threads) {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }
    
    crs_mt = (SPM_CRS_SYM_MT_TYPE *) malloc(ncpus * sizeof(spm_crs32_double_sym_mt_t));
    if (!crs_mt ) {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }

    cur_row = 0;
    elems_total = 0;
    for (i = 0; i < ncpus; i++) {
        elems_limit = (crs->nnz + crs->n - elems_total) / (ncpus - i);        
        
        spm_thread = spm_mt->spm_threads + i;
        spm_thread->cpu = cpus[i];
        spm_thread->spm = crs_mt + i;
        spm_thread->id = i;
        
        row_start = cur_row;
        crs_mt[i].row_start = cur_row;
		
        elems_total += elems_limit;
        while (crs->row_ptr[cur_row] + cur_row < elems_total)
            cur_row++;
        elems_total = cur_row + crs->row_ptr[cur_row];
        
        crs_mt[i].row_end = cur_row;
        crs_mt[i].nnz = elems_total;
        crs_mt[i].crs = crs;
    }
    
    assert(cur_row == crs->n);
    assert(elems_total == crs->nnz + crs->n);
    
    free(cpus);
    return spm_mt;
}

void SPM_CRS_SYM_MT_NAME(_destroy)(void *spm)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
    SPM_CRS_SYM_MT_TYPE *crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
    
    SPM_CRS_SYM_NAME(_destroy)(crs_mt->crs);
    free(crs_mt);
    free(spm_thread);
    free(spm_mt);
}

uint64_t SPM_CRS_SYM_MT_NAME(_size)(void *spm)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
    SPM_CRS_SYM_MT_TYPE *crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
    
    return SPM_CRS_SYM_NAME(_size)(crs_mt->crs);
}

void SPM_CRS_SYM_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
    SPM_CRS_SYM_MT_TYPE *crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm;
    ELEM_TYPE *x = in->elements;
    ELEM_TYPE *y = out->elements;
    ELEM_TYPE *dvalues = crs_mt->crs->dvalues;
    ELEM_TYPE *values = crs_mt->crs->values;
    SPM_CRS_SYM_IDX_TYPE *row_ptr = crs_mt->crs->row_ptr;
    SPM_CRS_SYM_IDX_TYPE *col_ind = crs_mt->crs->col_ind;
    const uint64_t row_start = crs_mt->row_start;
    const uint64_t row_end = crs_mt->row_end;
    register double yr;
    uint64_t i,j;
    
    ///> Parallel multiplications.
    for (i = row_start; i < row_end; i++) {
        yr = 0;
        for (j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            yr += values[j] * x[col_ind[j]];
            y[col_ind[j]] += values[j] * x[i];
        }
        yr += dvalues[i] * x[i];
        y[i] = yr;
    }
}

XSPMV_SYM_MT_METH_INIT(
    SPM_CRS_SYM_MT_NAME(_multiply),
    SPM_CRS_SYM_MT_NAME(_init_mmf),
    SPM_CRS_SYM_MT_NAME(_size),
    SPM_CRS_SYM_MT_NAME(_destroy),
    sizeof(ELEM_TYPE)
)

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
