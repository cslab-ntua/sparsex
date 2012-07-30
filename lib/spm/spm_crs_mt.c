/*
 * spm_crs_mt.c -- Multithreaded CRS implementation.
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <assert.h>
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "macros.h"
#include "mt_lib.h"
#include "spm_mt.h"
#include "spm_crs.h"
#include "spm_crs_mt.h"
#include "spmv_method.h"

// Set to true when in special benchmark context. It is used to implement
// the noxmiss benchmark (sequential access pattern in the x vector)
static int SpmBenchContext = 0;

static double arith_intensity(size_t nr_rows, size_t nr_nzeros)
{
    size_t elem_size = sizeof(ELEM_TYPE);
    size_t idx_size = sizeof(SPM_CRS_IDX_TYPE);
    return (2.*nr_nzeros) /
        (nr_nzeros*(elem_size + idx_size) + nr_rows*(idx_size + 2.*elem_size));
}

void *SPM_CRS_MT_NAME(_init_mmf_noxmiss)(char *mmf_file, uint64_t *rows_nr,
                                         uint64_t *cols_nr, uint64_t *nz_nr,
                                         void *metadata)
{
    SpmBenchContext = 1;
    return SPM_CRS_MT_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr, nz_nr,
                                      metadata);
}

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

	// Compute the matrix splits.
	size_t nnz_per_split = csr->nz / nr_threads;
	size_t curr_nnz = 0;
	size_t row_start = 0;
	size_t split_cnt = 0;
	int i;
	for (i = 0; i < nr_rows; i++) {
		curr_nnz += rowptr[i+1] - rowptr[i];
		if (curr_nnz >= nnz_per_split) {
			csr_mts[split_cnt].crs = csr;
			csr_mts[split_cnt].row_start = row_start;
			csr_mts[split_cnt].row_end = i + 1;
			ret->spm_threads[split_cnt].spm = &csr_mts[split_cnt];
			ret->spm_threads[split_cnt].spmv_fn =
			    SPM_CRS_MT_NAME(_multiply);
			ret->spm_threads[split_cnt].cpu = cpus[split_cnt];
			ret->spm_threads[split_cnt].row_start = row_start;
			ret->spm_threads[split_cnt].nr_rows = i + 1 - row_start;
			row_start = i + 1;
			curr_nnz = 0;
			++split_cnt;
		}
	}

	// Fill the last split.
	if (curr_nnz < nnz_per_split && split_cnt < nr_threads) {
		csr_mts[split_cnt].crs = csr;
		csr_mts[split_cnt].row_start = row_start;
		csr_mts[split_cnt].row_end = i;
		ret->spm_threads[split_cnt].spm = &csr_mts[split_cnt];
		ret->spm_threads[split_cnt].spmv_fn =
		    SPM_CRS_MT_NAME(_multiply);
		ret->spm_threads[split_cnt].cpu = cpus[split_cnt];
		ret->spm_threads[split_cnt].row_start = row_start;
		ret->spm_threads[split_cnt].nr_rows = i - row_start;
	}

	return ret;
}


void *SPM_CRS_MT_NAME(_init_mmf)(char *mmf_file, uint64_t *rows_nr,
                                 uint64_t *cols_nr, uint64_t *nz_nr,
                                 void *metadata)
{
	int i;
	unsigned int nr_cpus, *cpus;
	spm_mt_t *spm_mt;
	spm_mt_thread_t *spm_thread;
	SPM_CRS_MT_TYPE *crs_mt;
	unsigned long cur_row, elems_limit, elems_total=0;

	// Set affinity of the current thread.
	mt_get_options(&nr_cpus, &cpus);
	setaffinity_oncpu(cpus[0]);

	// Print the chosen cpus configuration given at command line.
	printf("MT_CONF:%d", cpus[0]);
	for (i = 1; i < nr_cpus; i++)
		printf(",%d", cpus[i]);
	printf("\n");

	SPM_CRS_TYPE *crs;
	if (SpmBenchContext)
		crs = SPM_CRS_NAME(_init_mmf_noxmiss)(mmf_file, rows_nr, cols_nr,
		                                      nz_nr, metadata);
	else
		crs = SPM_CRS_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr,
		                              nz_nr, metadata);

#if SPM_INTENSITY_MAP
	const char *map_file = "mat_intensity.out";
	FILE *fmap = fopen(map_file, "w");
	if (!fmap) {
		perror("WARNING: fopen() failed");
		fprintf(stderr, "File: %s\n", map_file);
	}
#endif

	spm_mt = malloc(sizeof(spm_mt_t));
	if (!spm_mt) {
		perror("malloc");
		exit(1);
	}

	spm_mt->nr_threads = nr_cpus;
	spm_mt->spm_threads = malloc(sizeof(spm_mt_thread_t)*nr_cpus);
	if (!spm_mt->spm_threads) {
		perror("malloc");
		exit(1);
	}

	crs_mt = malloc(sizeof(SPM_CRS_MT_TYPE)*nr_cpus);
	if (!crs_mt) {
		perror("malloc");
		exit(1);
	}

	double total_intensity = arith_intensity(crs->nrows, crs->nz);
	printf("Matrix flop:byte ratio: %.4f\n", total_intensity);
	for (i = 0, cur_row = 0; i < nr_cpus; i++) {
		double local_intensity __attribute__((unused)) = 0.0;
		elems_limit = (*nz_nr - elems_total) / (nr_cpus - i);
		spm_thread = spm_mt->spm_threads+ i;
		spm_thread->cpu = cpus[i];
		spm_thread->spm = crs_mt + i;

		unsigned long elems = 0;
		size_t sample_elems __attribute__ ((unused)) = 0;
		crs_mt[i].row_start = cur_row;
		while (1) {
			elems += crs->row_ptr[cur_row+1] - crs->row_ptr[cur_row];
#ifdef SPM_INTENSITY_MAP
			sample_elems += crs->row_ptr[cur_row+1] - crs->row_ptr[cur_row];
			if (cur_row && cur_row % 1000 == 0) {
				local_intensity = arith_intensity(1000, sample_elems);
				if (fmap)
					fprintf(fmap, "%zd %lf\n", cur_row, local_intensity);
				sample_elems = 0;
			}

#endif
			cur_row++;
//			  if (cur_row % row_limit == 0 || cur_row == *rows_nr)
			if (elems >= elems_limit)
				break;
		}

		elems_total += elems;
		crs_mt[i].row_end = cur_row;
		crs_mt[i].crs = crs;
		int nr_rows_part = crs_mt[i].row_end - crs_mt[i].row_start;
		int nr_nzeros_part = elems;
		printf("Partition info (id, nr_rows, nr_nzeros, flop:byte): (%d,%d,%d,%.4f)\n",
		       i, nr_rows_part, nr_nzeros_part,
		       arith_intensity(nr_rows_part, nr_nzeros_part));
	}

#ifdef SPM_INTENSITY_MAP
    if (fmap)
        fclose(fmap);
#endif

	free(cpus);
	return spm_mt;
}

void SPM_CRS_MT_NAME(_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *) spm_thread->spm;

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

	for (i = row_start; i < row_end; i++) {
		yr = (ELEM_TYPE) 0;
		for (j = row_ptr[i]; j < row_ptr[i+1]; j++)
			yr += (values[j] * x[col_ind[j]]);

		y[i] = yr;
	}
}

#ifdef SPM_BENCH
/*
 *  Benchmark routines
 *
 *  All benchmarks assume the special noxmiss initialization of the matrix.
 *  See `SpmBenchContext' in `spm_crs.c' file.
 */
void SPM_CRS_MT_NAME(_noxmiss_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_NAME(_multiply)(spm, in, out);
}

void SPM_CRS_MT_NAME(_nocolind_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	SPM_CRS_IDX_TYPE *row_ptr = crs_mt->crs->row_ptr;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	register ELEM_TYPE yr;
	unsigned long i,j, k;

	for (i = row_start; i < row_end; i++) {
		yr = (ELEM_TYPE) 0;
		for (j = row_ptr[i], k = 0; j < row_ptr[i+1]; j++, k++)
			yr += (values[j] * x[k]);

		y[i] = yr;
	}
}

void SPM_CRS_MT_NAME(_norowptr_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	SPM_CRS_IDX_TYPE *col_ind = crs_mt->crs->col_ind;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	register ELEM_TYPE yr;
	unsigned long i,j;

	size_t nnz = crs_mt->crs->nz;
	size_t nrows = crs_mt->crs->nrows;
	size_t row_size = nnz / nrows + (nnz % nrows != 0);
	for (i = row_start; i < row_end; i++) {
		size_t rs = i*row_size;
		yr = (ELEM_TYPE) 0;
		for (j = 0; j < row_size; j++)
			yr += (values[rs+j] * x[col_ind[rs+j]]);

		y[i] = yr;
	}
}

void SPM_CRS_MT_NAME(_norowptr_nocolind_multiply)(void *spm, VECTOR_TYPE *in, VECTOR_TYPE *out)
{
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *)spm;
	ELEM_TYPE *x = in->elements;
	ELEM_TYPE *y = out->elements;
	ELEM_TYPE *values = crs_mt->crs->values;
	const unsigned long row_start = crs_mt->row_start;
	const unsigned long row_end = crs_mt->row_end;
	register ELEM_TYPE yr;
	unsigned long i,j,k;

	size_t nnz = crs_mt->crs->nz;
	size_t nrows = crs_mt->crs->nrows;
	size_t row_size = nnz / nrows + (nnz % nrows != 0);
	for (i = row_start; i < row_end; i++) {
		size_t rs = i*row_size;
		yr = (ELEM_TYPE) 0;
		for (j = 0, k = 0; j < row_size; j++, k++)
			yr += values[rs+j]*x[k];

		y[i] = yr;
	}
}

XSPMV_MT_METH_INIT(
	SPM_CRS_MT_NAME(_noxmiss_multiply),
	SPM_CRS_MT_NAME(_init_mmf_noxmiss),
	SPM_CRS_MT_NAME(_size),
	SPM_CRS_MT_NAME(_destroy),
	sizeof(ELEM_TYPE)
)

XSPMV_MT_METH_INIT(
	SPM_CRS_MT_NAME(_norowptr_multiply),
	SPM_CRS_MT_NAME(_init_mmf_noxmiss),
	SPM_CRS_MT_NAME(_size),
	SPM_CRS_MT_NAME(_destroy),
	sizeof(ELEM_TYPE)
)

XSPMV_MT_METH_INIT(
	SPM_CRS_MT_NAME(_nocolind_multiply),
	SPM_CRS_MT_NAME(_init_mmf_noxmiss),
	SPM_CRS_MT_NAME(_size),
	SPM_CRS_MT_NAME(_destroy),
	sizeof(ELEM_TYPE)
)

XSPMV_MT_METH_INIT(
	SPM_CRS_MT_NAME(_norowptr_nocolind_multiply),
	SPM_CRS_MT_NAME(_init_mmf_noxmiss),
	SPM_CRS_MT_NAME(_size),
	SPM_CRS_MT_NAME(_destroy),
	sizeof(ELEM_TYPE)
)

#endif  /* SPM_BENCH */

// Multiplication for one-based indexing of rows and column indices.
void SPM_CRS_MT_NAME(_multiply_base_one)(void *spm, VECTOR_TYPE *in,
                                         VECTOR_TYPE *out)
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

	for (i = row_start; i < row_end; i++) {
		yr = (ELEM_TYPE) 0;
		for (j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
			yr += (values[j] * x[col_ind[j]-1]);

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

void *SPM_CRS_MT_NAME(_numa_init_mmf)(char *mmf_file, uint64_t *rows_nr,
                                      uint64_t *cols_nr, uint64_t *nz_nr,
                                      void *metadata)
{
	spm_mt_t *spm_mt = SPM_CRS_MT_NAME(_init_mmf)(mmf_file, rows_nr, cols_nr,
	                                              nz_nr, metadata);

	int nr_threads = spm_mt->nr_threads;
	size_t *values_parts = malloc(nr_threads*sizeof(*values_parts));
	size_t *rowptr_parts = malloc(nr_threads*sizeof(*rowptr_parts));
	size_t *colind_parts = malloc(nr_threads*sizeof(*colind_parts));
	int *nodes = malloc(nr_threads*sizeof(*nodes));

	// Reallocate data structures in a numa-aware fashion.
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

	// Sanity check.
	assert(crs);

	SPM_CRS_IDX_TYPE *new_rowptr =
	    alloc_interleaved((crs->nrows+1)*sizeof(*crs->row_ptr), rowptr_parts,
	                      nr_threads, nodes);

	SPM_CRS_IDX_TYPE *new_colind =
	    alloc_interleaved(crs->nz*sizeof(*crs->col_ind), colind_parts,
	                      nr_threads, nodes);

	ELEM_TYPE *new_values =
	    alloc_interleaved(crs->nz*sizeof(*crs->values), values_parts,
	                      nr_threads, nodes);

	// Copy the old data to the new one.
	memcpy(new_rowptr, crs->row_ptr, (crs->nrows+1)*sizeof(*crs->row_ptr));
	memcpy(new_colind, crs->col_ind, crs->nz*sizeof(*crs->col_ind));
	memcpy(new_values, crs->values, crs->nz*sizeof(*crs->values));

	// Free old data and replace with the new one.
	free(crs->row_ptr);
	free(crs->col_ind);
	free(crs->values);

	crs->row_ptr = new_rowptr;
	crs->col_ind = new_colind;
	crs->values = new_values;

	// Check for allocation errors.
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

	// SPM_CRS_MT_NAME(_numa_make_col_map)(spm_mt);

	// Free the auxiliaries.
	free(rowptr_parts);
	free(colind_parts);
	free(values_parts);
	free(nodes);

	return spm_mt;
}

/* void SPM_CRS_MT_NAME(_numa_make_col_map)(void *spm)
{
	int i;
	uint64_t j, r, c;
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *) spm_thread->spm;
	SPM_CRS_TYPE *crs = crs_mt->crs;
	uint64_t *map;

	for (i = 0; i < spm_mt->nr_threads; i++) {
		map = (uint64_t *) malloc(crs->ncols*sizeof(uint64_t));
		for (c = 0; c < crs->ncols; c++)
			map[c] = 0;

		for (r = crs_mt[i].row_start; r < crs_mt[i].row_end; r++) {
			for (j = crs->row_ptr[r]; j < crs->row_ptr[r+1]; j++) {
				c = crs->col_ind[j];
				map[c]++;
			}
		}
		spm_thread[i].col_map = map;
	}
}*/

void SPM_CRS_MT_NAME(_numa_destroy)(void *spm)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads;
	SPM_CRS_MT_TYPE *crs_mt = (SPM_CRS_MT_TYPE *) spm_thread->spm;
	SPM_CRS_TYPE *crs = crs_mt->crs;

	free_interleaved(crs->row_ptr, (crs->nrows+1)*sizeof(*crs->row_ptr));
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

void SPM_CRS_MT_NAME(_numa_multiply)(void *spm, VECTOR_TYPE *in,
                                     VECTOR_TYPE *out)
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
