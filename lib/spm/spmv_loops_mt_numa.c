/*
 * spmv_loops_mt_numa.c -- NUMA-aware SpMV implementation
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <assert.h>
#include <stdio.h>
#include <pthread.h>
#include <numa.h>
#include <numaif.h>
#include "spm_mt.h"
#include "mt_lib.h"
#include "spmv_loops_mt_numa.h"
#include "tsc.h"
#include "timer.h"
#include "vector.h"

static VECTOR_TYPE *y = NULL;
static VECTOR_TYPE *x = NULL;
static pthread_barrier_t barrier;
static unsigned long loops_nr = 0;
static float secs = 0.0;

static void *do_spmv_thread_main(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
	SPMV_NAME(_fn_t) *spmv_mt_fn = spm_mt_thread->spmv_fn;
	setaffinity_oncpu(spm_mt_thread->cpu);

	int i;
	tsc_t total_tsc, thread_tsc;

	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
	for (i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
		tsc_start(&thread_tsc);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
	}

	tsc_pause(&total_tsc);
	spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
	secs = tsc_getsecs(&total_tsc);
	tsc_shut(&thread_tsc);
	tsc_shut(&total_tsc);

	return (void *) 0;
}

static void *do_spmv_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
	SPMV_NAME(_fn_t) *spmv_mt_fn = spm_mt_thread->spmv_fn;
	setaffinity_oncpu(spm_mt_thread->cpu);

	int i;
	tsc_t thread_tsc;

	tsc_init(&thread_tsc);
	for (i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
		tsc_start(&thread_tsc);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
	}

	spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
	tsc_shut(&thread_tsc);

	return NULL;
}

float SPMV_NAME(_bench_mt_loop_numa)(spm_mt_t *spm_mt, unsigned long loops,
                                     unsigned long rows_nr,
                                     unsigned long cols_nr,
                                     SPMV_NAME(_fn_t) *fn)
{
	int err, i;
	pthread_t *tids;

	setaffinity_oncpu(spm_mt->spm_threads[0].cpu);
	loops_nr = loops;

	err = pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads);
	if (err){
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = malloc(sizeof(pthread_t)*spm_mt->nr_threads);
	if (!tids){
		perror("malloc");
		exit(1);
	}

	/*
	int pagesize = numa_pagesize();
	uint64_t xparts_size = (cols_nr * sizeof(ELEM_TYPE)) / pagesize;
    
	if ((cols_nr * sizeof(ELEM_TYPE)) % pagesize != 0)
		xparts_size++;

	size_t *xparts = malloc(sizeof(*xparts)*xparts_size);
	int *xnodes = malloc(sizeof(*xnodes)*xparts_size);
	if (!xparts || !xnodes) {
		perror("malloc");
		exit(1);
	}

	int nr_nodes = numa_num_configured_nodes();
	uint64_t *count = (uint64_t *) malloc(nr_nodes*sizeof(uint64_t));
	uint64_t xpart, k, l, max, chosen_node;

	for (xpart = 0; xpart < xparts_size; xpart++) {
		for (k = 0; k < nr_nodes; k++)
			count[k] = 0;

		for (i = 0; i < spm_mt->nr_threads; i++) {
			spm_mt_thread_t *spm = spm_mt->spm_threads + i;
			int node = spm->node;

			for (l = xpart*pagesize/sizeof(ELEM_TYPE);
			     l < (xpart+1)*pagesize/sizeof(ELEM_TYPE) && l < cols_nr; l++)
				count[node] += spm->col_map[l];
		}

		chosen_node = 0;
		max = 0;
		for (k = 0; k < nr_nodes; k++) {
			if (count[k] > max) {
				max = count[k];
				chosen_node = k;
			}
		}

		if (xpart != xparts_size - 1)
			xparts[xpart] = pagesize;
		else
			xparts[xpart] = cols_nr*sizeof(ELEM_TYPE) -
			                (xparts_size - 1) * pagesize;
		xnodes[xpart] = chosen_node;
	}
	*/
	size_t *parts = malloc(sizeof(*parts)*spm_mt->nr_threads);
	int *nodes = malloc(sizeof(*nodes)*spm_mt->nr_threads);
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	for (i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;

		parts[i] = spm->nr_rows * sizeof(ELEM_TYPE);
		nodes[i] = spm->node;
		if (fn)
			spm->spmv_fn = fn;
    }

	// Allocate an interleaved x.
	int alloc_err = 0;

	x = VECTOR_NAME(_create_interleaved)(rows_nr, parts, spm_mt->nr_threads,
	                                     nodes);
	VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);
	alloc_err = check_interleaved(x->elements, parts, spm_mt->nr_threads,
	                              nodes);
	print_alloc_status("input vector", alloc_err);

	// Allocate an interleaved y.
	y = VECTOR_NAME(_create_interleaved)(rows_nr, parts, spm_mt->nr_threads,
	                                     nodes);
	VECTOR_NAME(_init)(y, 0);
	alloc_err = check_interleaved(y->elements, parts, spm_mt->nr_threads,
	                              nodes);
	print_alloc_status("output vector", alloc_err);

	// Run benchmark.
	pthread_create(tids, NULL, do_spmv_thread_main, spm_mt->spm_threads);
	for (i = 1; i < spm_mt->nr_threads; i++)
		pthread_create(tids + i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (i = 0; i < spm_mt->nr_threads; i++)
		pthread_join(tids[i], NULL);

	// Destroy vectors.
	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	pthread_barrier_destroy(&barrier);
	// free(count);
	// free(xparts);
	// free(xnodes);
	free(parts);
	free(nodes);
	free(tids);

	return secs;
}

void SPMV_NAME(_check_mt_loop_numa)(void *spm_serial, spm_mt_t *spm_mt,
                                    SPMV_NAME(_fn_t) *fn, unsigned long loops,
                                    unsigned long rows_nr,
                                    unsigned long cols_nr,
                                    SPMV_NAME(_fn_t) *mt_fn)
{
	int err, i;
	pthread_t *tids;
	VECTOR_TYPE *y2;

	loops_nr = loops;
	err = pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads + 1);
	if (err) {
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = malloc(sizeof(pthread_t)*spm_mt->nr_threads);
	if (!tids) {
		perror("malloc");
		exit(1);
	}

	size_t *parts = malloc(sizeof(*parts)*spm_mt->nr_threads);
	int *nodes = malloc(sizeof(*nodes)*spm_mt->nr_threads);
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	for (i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);

		parts[i] = spm->nr_rows * sizeof(ELEM_TYPE);
		nodes[i] = spm->node;
		if (mt_fn)
			spm->spmv_fn = mt_fn;
	}

	int node = spm_mt->spm_threads[0].node;

	x = VECTOR_NAME(_create_onnode)(cols_nr, node);
	VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);

	y = VECTOR_NAME(_create_interleaved)(rows_nr, parts, spm_mt->nr_threads,
	                                     nodes);
	y2 = VECTOR_NAME(_create)(rows_nr);
	VECTOR_NAME(_init)(y, 21);
	VECTOR_NAME(_init)(y2, 0);

	for (i = 0; i < spm_mt->nr_threads; i++)
		pthread_create(tids + i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		// Use the prototype of x.
		fn(spm_serial, x, y2);
		if (VECTOR_NAME(_compare)(y2, y) < 0)
			exit(1);
	}

	for (i = 0; i < spm_mt->nr_threads; i++)
		pthread_join(tids[i], NULL);

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	VECTOR_NAME(_destroy)(y2);

	pthread_barrier_destroy(&barrier);
	free(parts);
	free(nodes);
	free(tids);
}
