/*
 * spm_crs_mt.c -- Multithreaded CRS implementation.
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C)      2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdint.h>
#include <assert.h>

#include "macros.h"
#include "mt_lib.h"
#include "spm_mt.h"
#include "spm_crs.h"
#include "spm_crs_mt.h"
#include "spmv_method.h"

spm_mt_t *SPM_CRS_MT_NAME(_get_spm)(SPM_CRS_IDX_TYPE *rowptr,
                                    SPM_CRS_IDX_TYPE *colind,
                                    ELEM_TYPE *values,
                                    SPM_CRS_IDX_TYPE nr_rows,
                                    unsigned int nr_threads,
                                    unsigned int *cpus)
{
    SPM_CRS_MT_TYPE *csr_mts = malloc(nr_threads*sizeof(*csr_mts));
    SPM_CRS_TYPE *csr = malloc(sizeof(*csr));
    csr->row_ptr = rowptr;
    csr->col_ind = colind;
    csr->values = values;
    csr->nz = rowptr[nr_rows] - 1;
    csr->nrows = csr->ncols = nr_rows;

    spm_mt_t *ret = malloc(sizeof(*ret));
    ret->nr_threads = nr_threads;
    ret->spm_threads = malloc(nr_threads*sizeof(*ret->spm_threads));

    // Compute the matrix splits
    size_t nnz_per_split = csr->nz / nr_threads;
    size_t curr_nnz = 0;
    size_t row_start = 0;
    size_t split_cnt = 0;
    int i;
    for (i = 0; i < nr_rows; ++i) {
        curr_nnz += rowptr[i+1] - rowptr[i];
        if (curr_nnz >= nnz_per_split) {
            csr_mts[split_cnt].crs = csr;
            csr_mts[split_cnt].row_start = row_start;
            csr_mts[split_cnt].row_end = i+1;
            ret->spm_threads[split_cnt].spm = &csr_mts[split_cnt];
            ret->spm_threads[split_cnt].spmv_fn =
                SPM_CRS_MT_NAME(_multiply_base_one);
            ret->spm_threads[split_cnt].cpu = cpus[split_cnt];
            ret->spm_threads[split_cnt].row_start = row_start;
            ret->spm_threads[split_cnt].nr_rows = i + 1 - row_start;
            row_start = i + 1;
            curr_nnz = 0;
            ++split_cnt;
        }
    }

    // fill the last split
    if (curr_nnz < nnz_per_split && split_cnt < nr_threads) {
        csr_mts[split_cnt].crs = csr;
        csr_mts[split_cnt].row_start = row_start;
        csr_mts[split_cnt].row_end = i;
        ret->spm_threads[split_cnt].spm = &csr_mts[split_cnt];
        ret->spm_threads[split_cnt].spmv_fn =
            SPM_CRS_MT_NAME(_multiply_base_one);
        ret->spm_threads[split_cnt].cpu = cpus[split_cnt];
        ret->spm_threads[split_cnt].row_start = row_start;
        ret->spm_threads[split_cnt].nr_rows = i - row_start;
    }

    return ret;
}


void *SPM_CRS_MT_NAME(_init_mmf)(char *mmf_file,
                                 uint64_t *rows_nr, uint64_t *cols_nr,
                                 uint64_t *nz_nr, void *metadata)
{
	int i;
	unsigned int nr_cpus, *cpus;
	spm_mt_t *spm_mt;
	spm_mt_thread_t *spm_thread;
	SPM_CRS_MT_TYPE *crs_mt;
	unsigned long cur_row, elems_limit, elems_total=0;

	// set affinity of the current thread
	mt_get_options(&nr_cpus, &cpus);
	setaffinity_oncpu(cpus[0]);

	printf("MT_CONF:%d", cpus[0]);
	for (i = 1; i < nr_cpus; i++)
		printf(",%d", cpus[i]);
	printf("\n");

	SPM_CRS_TYPE *crs;
	crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr, metadata);

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

	crs_mt = malloc(sizeof(SPM_CRS_MT_TYPE)*nr_cpus);
	if ( !crs_mt ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	for (i=0, cur_row=0; i<nr_cpus; i++){
		elems_limit = (*nz_nr - elems_total) / (nr_cpus - i);
		spm_thread = spm_mt->spm_threads+ i;
		spm_thread->cpu = cpus[i];
		spm_thread->spm = crs_mt + i;

		unsigned long elems;
		crs_mt[i].row_start = cur_row;
		for (elems=0; ; ){
			elems += crs->row_ptr[cur_row+1] - crs->row_ptr[cur_row];
			cur_row++;
			//printf("i:%d nr_cpus:%d cur_row:%lu rows_nr:%lu elems:%lu elems_limit:%lu\n", i, nr_cpus, cur_row, *rows_nr, elems, elems_limit);
#if 0
			if (i != (nr_cpus -1)){
				if ( elems >= elems_limit )
					break;
			} else {
				if (cur_row == *rows_nr)
					break;
			}
#endif
			if (elems >= elems_limit){
				break;
			}
		}

		elems_total += elems;
		crs_mt[i].row_end = cur_row;
		crs_mt[i].crs = crs;
	}

	free(cpus);
	return spm_mt;
}

void SPM_CRS_MT_NAME(_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *)spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm_thread->spm;
	SPM_CRS_NAME(_destroy)(crs_mt->crs);
	free(crs_mt);
	free(spm_thread);
	free(spm_mt);
}

uint64_t SPM_CRS_MT_NAME(_size)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *) spm_thread->spm;
	return SPM_CRS_NAME(_size)(crs_mt->crs);
}

void SPM_CRS_MT_NAME(_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	SPM_CRS_IDX_TYPE *row_ptr = crs_mt->crs->row_ptr;
	SPM_CRS_IDX_TYPE *col_ind = crs_mt->crs->col_ind;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	register ELEM_TYPE yr;

	unsigned long i,j;
	for (i=row_start; i<row_end; i++){
		yr = (ELEM_TYPE)0;
		for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
			yr += (values[j] * x[col_ind[j]]);
		}
		y[i] = yr;
	}
}

/*
 *  Multiplication for one-based indexing of rows and column indices
 */
void SPM_CRS_MT_NAME(_multiply_base_one)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	SPM_CRS_IDX_TYPE *row_ptr = crs_mt->crs->row_ptr;
	SPM_CRS_IDX_TYPE *col_ind = crs_mt->crs->col_ind;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	register ELEM_TYPE yr;

	unsigned long i,j;
	for (i=row_start; i<row_end; i++){
		yr = (ELEM_TYPE)0;
		for (j=row_ptr[i]-1; j<row_ptr[i+1]-1; j++) {
			yr += (values[j] * x[col_ind[j]-1]);
		}
		y[i] = yr;
	}
}

XSPMV_MT_METH_INIT(
 SPM_CRS_MT_NAME(_multiply),
 SPM_CRS_MT_NAME(_init_mmf),
 SPM_CRS_MT_NAME(_size),
 SPM_CRS_MT_NAME(_destroy),
 sizeof(ELEM_TYPE)
)

#ifdef SPM_NUMA

#include <numa.h>
#include "numa_util.h"

void *SPM_CRS_MT_NAME(_numa_init_mmf)(char *mmf_file,
                                      uint64_t *rows_nr, uint64_t *cols_nr,
                                      uint64_t *nz_nr, void *metadata)
{
	spm_mt_t *spm_mt = SPM_CRS_MT_NAME(_init_mmf)(mmf_file,
	                                              rows_nr, cols_nr,
	                                              nz_nr, metadata);

	int nr_threads = spm_mt->nr_threads;
	size_t *values_parts = malloc(nr_threads*sizeof(*values_parts));
	size_t *rowptr_parts = malloc(nr_threads*sizeof(*rowptr_parts));
	size_t *colind_parts = malloc(nr_threads*sizeof(*colind_parts));
	int *nodes = malloc(nr_threads*sizeof(*nodes));

	// just reallocate in a numa-aware fashion the data structures
	int i;
	SPM_CRS_TYPE *crs = NULL;
	for (i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
		SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *) spm_thread->spm;
		crs = crs_mt->crs;
		SPM_CRS_IDX_TYPE row_start = crs_mt->row_start;
		SPM_CRS_IDX_TYPE row_end = crs_mt->row_end;
		SPM_CRS_IDX_TYPE vstart = crs->row_ptr[crs_mt->row_start];
		SPM_CRS_IDX_TYPE vend = crs->row_ptr[crs_mt->row_end];
		rowptr_parts[i] = (row_end-row_start)*sizeof(*crs->row_ptr);
		colind_parts[i] = (vend-vstart)*sizeof(*crs->col_ind);
		values_parts[i] = (vend-vstart)*sizeof(*crs->values);
		nodes[i] = numa_node_of_cpu(spm_thread->cpu);
		spm_thread->node = nodes[i];
		spm_thread->row_start = row_start;
		spm_thread->nr_rows = row_end - row_start;
	}
	rowptr_parts[nr_threads-1] += sizeof(*crs->row_ptr); 

	// sanity check (+ get rid of compiler warning about uninit. variable)
	assert(crs);

	SPM_CRS_IDX_TYPE *new_rowptr = alloc_interleaved((crs->nrows+1)*sizeof(*crs->row_ptr),
	                                                 rowptr_parts, nr_threads,
	                                                 nodes);
	SPM_CRS_IDX_TYPE *new_colind = alloc_interleaved(crs->nz*sizeof(*crs->col_ind),
	                                                 colind_parts, nr_threads,
	                                                 nodes);
	ELEM_TYPE *new_values = alloc_interleaved(crs->nz*sizeof(*crs->values),
	                                          values_parts, nr_threads,
	                                          nodes);

	// copy old data to the new one
	memcpy(new_rowptr, crs->row_ptr, (crs->nrows+1)*sizeof(*crs->row_ptr));
	memcpy(new_colind, crs->col_ind, crs->nz*sizeof(*crs->col_ind));
	memcpy(new_values, crs->values, crs->nz*sizeof(*crs->values));

	// free old data and replace with the new one
	free(crs->row_ptr);
	free(crs->col_ind);
	free(crs->values);
	crs->row_ptr = new_rowptr;
	crs->col_ind = new_colind;
	crs->values = new_values;

	int alloc_err;
	alloc_err = check_interleaved((void *) crs->row_ptr, rowptr_parts,
	                              nr_threads, nodes);
	print_alloc_status("CSR rowptr", alloc_err);
	alloc_err = check_interleaved((void *) crs->col_ind,  colind_parts,
	                              nr_threads, nodes);
	print_alloc_status("CSR colind", alloc_err);
	alloc_err = check_interleaved((void *) crs->values, values_parts,
	                              nr_threads, nodes);
	print_alloc_status("CSR values", alloc_err);

	// free the auxiliaries
	free(rowptr_parts);
	free(colind_parts);
	free(values_parts);
	free(nodes);
	return spm_mt;
}

void SPM_CRS_MT_NAME(_numa_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *) spm_thread->spm;
	SPM_CRS_TYPE *crs = crs_mt->crs;

	free_interleaved(crs->col_ind, crs->nz*sizeof(*crs->col_ind));
	free_interleaved(crs->values, crs->nz*sizeof(*crs->values));
	free(crs);
	free(crs_mt);
	free(spm_thread);
	free(spm_mt);
}

uint64_t SPM_CRS_MT_NAME(_numa_size)(void *spm)
{
	return SPM_CRS_MT_NAME(_size)(spm);
}

void SPM_CRS_MT_NAME(_numa_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_NAME(_multiply)(spm, in, out);
}

XSPMV_MT_METH_INIT(
 SPM_CRS_MT_NAME(_numa_multiply),
 SPM_CRS_MT_NAME(_numa_init_mmf),
 SPM_CRS_MT_NAME(_numa_size),
 SPM_CRS_MT_NAME(_numa_destroy),
 sizeof(ELEM_TYPE)
)

#endif /* SPM_NUMA */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
