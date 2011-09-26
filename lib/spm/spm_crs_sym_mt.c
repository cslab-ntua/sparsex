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
    
    crs_mt = (SPM_CRS_SYM_MT_TYPE *) malloc(ncpus * sizeof(SPM_CRS_SYM_MT_TYPE));
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
    
    SPM_CRS_SYM_MT_NAME(_make_map)(spm_mt);
    
    assert(cur_row == crs->n);
    assert(elems_total == crs->nnz + crs->n);
    
    free(cpus);
    return spm_mt;
}

void SPM_CRS_SYM_MT_NAME(_make_map)(void *spm) {
    unsigned int i;
    uint64_t j, k;
    uint64_t start, end;
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    unsigned int ncpus = spm_mt->nr_threads;
    spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
    SPM_CRS_SYM_MT_TYPE *crs_mt;
    SPM_CRS_SYM_TYPE *crs;
    SPM_CRS_SYM_IDX_TYPE *row_ptr;
    SPM_CRS_SYM_IDX_TYPE *col_ind;
    map_t *final_map;
    unsigned int total_count = 0;
    unsigned int temp_count, limit;
    int **initial_map;
    unsigned int *count;
    
    ///> Init parameters.
    spm_thread = spm_mt->spm_threads;
    crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
    crs = crs_mt->crs;
    row_ptr = crs->row_ptr;
    col_ind = crs->col_ind;
    
    initial_map = (int **) malloc(ncpus * sizeof(int *));
    for (i = 0; i < ncpus; i++)
        initial_map[i] = (int *) malloc(crs->n * sizeof(int));
    count = (unsigned int *) malloc(crs->n * sizeof(unsigned int));
    
    for (i = 0; i < ncpus; i++)
        for (j = 0; j < crs->n; j++)
            initial_map[i][j] = 0;
    for (j = 0; j < crs->n; j++)
        count[j] = 0;
        
    ///> Create initial map.  
    for (i = 0; i < ncpus; i++) {
        spm_thread = spm_mt->spm_threads + i;
        crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
        start = crs_mt->row_start;
        end = crs_mt->row_end;
        
        for (j = start; j < end; j++) {
            for (k = row_ptr[j]; col_ind[k] < start && k < row_ptr[j+1]; k++) {
                if (initial_map[i][col_ind[k]] == 0) {
                    initial_map[i][col_ind[k]] = 1;
                    count[col_ind[k]]++;
                }
            }
        }
    }
    
    ///> Print initial map.
    /*
    for (i = 0; i < ncpus; i++) {
        for (j = 0; j < crs->n; j++)
            printf("%d ", initial_map[i][j]);
        printf("\n");
    }
    
    for (j = 0; j < crs->n; j++)
        printf("%u ", count[j]);
    printf("\n");
    
    temp_count = 0;
    for (i = 0; i < ncpus; i++) {
        spm_thread = spm_mt->spm_threads + i;
        crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
        start = crs_mt->row_start;
        for (j = 0; j < start; j++)
            if (initial_map[i][j] == 1)
                temp_count++;
        total_count += start;
    }
    printf("%lf\n", (double) temp_count / (double) total_count);
    total_count = 0;
    */
    
    ///> Create final map.
    start = 0;
    end = 0;
    for (j = 0; j < crs->n; j++)
        total_count += count[j];
        
    for (i = 0; i < ncpus - 1; i++) {
        spm_thread = spm_mt->spm_threads + i;
        crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
        final_map = (map_t *) malloc(sizeof(map_t));
        
        limit = (total_count) / (ncpus - i);
        temp_count = 0;
        while (temp_count < limit) {
            temp_count += count[end];
            end++;
        }
        
        final_map->length = temp_count;
        final_map->cpus = (uint32_t *) malloc(temp_count * sizeof(uint32_t));
        final_map->elems_pos = 
            (uint32_t *) malloc(temp_count * sizeof(uint32_t));
                
        temp_count = 0;
        for (k = start; k < end; k++) {
            for (j = 0; j < ncpus; j++) {
                if (initial_map[j][k] == 1) {
                    final_map->cpus[temp_count] = j;
                    final_map->elems_pos[temp_count++] = k;
                }
            }
        }
        spm_thread->map = final_map;
        
        start = end;
        total_count -= temp_count;
    }
    end = crs->n;
    
    spm_thread = spm_mt->spm_threads + ncpus - 1;
    crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
    final_map = (map_t *) malloc(sizeof(map_t));
    
    temp_count = 0;
    for (j = start; j < end; j++)
        temp_count += count[j];
    total_count -= temp_count;
    assert(total_count == 0);
    final_map->length = temp_count;
    final_map->cpus = (uint32_t *) malloc(temp_count * sizeof(uint32_t));
    final_map->elems_pos = 
        (uint32_t *) malloc(temp_count * sizeof(uint32_t));
    
    temp_count = 0;
    for (k = start; k < end; k++) {
        for (j = 0; j < ncpus; j++) {
            if (initial_map[j][k] == 1) {
                final_map->cpus[temp_count] = j;
                final_map->elems_pos[temp_count++] = k;
            }
        }
    }
    spm_thread->map = final_map;
    
    ///> Print map.
    
    /*
    for (i = 0; i < ncpus; i++) {
        spm_thread = spm_mt->spm_threads + i;
        final_map = spm_thread->map;
        
        printf("Thread %d\n", i);
        for (j = 0; j < final_map->length; j++)
            printf("(%u, %u)\n", final_map->cpus[j], final_map->elems_pos[j]);
    }
    */
    
    printf("Map Size %lu\n", SPM_CRS_SYM_MT_NAME(_map_size)(spm_mt));
    
    ///> Free parameters.
    for (i = 0; i < ncpus; i++)
        free(initial_map[i]);
    free(initial_map);
    free(count);
}

void SPM_CRS_SYM_MT_NAME(_destroy)(void *spm)
{
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
    SPM_CRS_SYM_MT_TYPE *crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm_thread->spm;
    
    SPM_CRS_SYM_NAME(_destroy)(crs_mt->crs);
    free(spm_thread->map->cpus);
    free(spm_thread->map->elems_pos);
    free(spm_thread->map);
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

uint64_t SPM_CRS_SYM_MT_NAME(_map_size)(void *spm)
{
    unsigned int i;
    spm_mt_t *spm_mt = (spm_mt_t *) spm;
    spm_mt_thread_t *spm_thread;
    uint64_t size = 0;
    
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        size += spm_thread->map->length * sizeof(uint32_t);
        size += spm_thread->map->length * sizeof(UINT_TYPE(SPM_CRS_BITS));
        size += spm_thread->map->length * sizeof(ELEM_TYPE);
    }
    return size;
}

void SPM_CRS_SYM_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, 
                                    VECTOR_TYPE *out, VECTOR_TYPE *temp)
{
    SPM_CRS_SYM_MT_TYPE *crs_mt = (SPM_CRS_SYM_MT_TYPE *) spm;
    ELEM_TYPE *x = in->elements;
    ELEM_TYPE *y = out->elements;
    ELEM_TYPE *t = temp->elements;
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
        yr = (ELEM_TYPE) 0;
        for (j = row_ptr[i]; col_ind[j] < row_start && j < row_ptr[i+1]; j++) {
            yr += values[j] * x[col_ind[j]];
            t[col_ind[j]] += values[j] * x[i];
        }
        for ( ; j < row_ptr[i+1]; j++) {
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
