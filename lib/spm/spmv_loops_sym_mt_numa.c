/*
 * spmv_loops_mt_numa.c -- NUMA-aware SpMV implementation for symmetric matrices
 *
 * Copyright (C) 2010-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011,      Theodoros Gkountouvas
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
#include "spmv_loops_sym_mt_numa.h"
#include "tsc.h"
#include "vector.h"

static VECTOR_TYPE *y = NULL;
static VECTOR_TYPE **temp = NULL;
static pthread_barrier_t barrier;
static unsigned long nloops = 0;
static unsigned long ncpus = 0;
static unsigned long n = 0;
static float secs = 0.0;

static void *do_spmv_thread(void *arg)
{
    
    spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
    SPMV_NAME(_sym_fn_t) *spmv_mt_sym_fn = spm_mt_thread->spmv_fn;
    int id = spm_mt_thread->id;
    int i; // j, start, end;
    
    setaffinity_oncpu(spm_mt_thread->cpu);
    // start = (id * n) / ncpus;
    // end = ((id + 1) * n) / ncpus;
    
    for (i = 0; i < nloops; i++) {
        /*
	if (id != 0)
            VECTOR_NAME(_init)(temp[id], 0);
        */
	VECTOR_NAME(_init_from_map)(temp, 0, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
        spmv_mt_sym_fn(spm_mt_thread->spm, spm_mt_thread->data, y, temp[id]);
        pthread_barrier_wait(&barrier);
        /*
	for (j = 1; j < ncpus; j++)
            VECTOR_NAME(_add_part)(y, temp[j], y, start, end);
        */
	VECTOR_NAME(_add_from_map)(y, temp, y, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
    }

    return NULL;
}

static void *do_spmv_thread_main(void *arg)
{
    spm_mt_thread_t *spm_mt_thread = arg;
    SPMV_NAME(_sym_fn_t) *spmv_mt_sym_fn = spm_mt_thread->spmv_fn;
    int id = spm_mt_thread->id;
    
    setaffinity_oncpu(spm_mt_thread->cpu);

    tsc_t tsc;
    tsc_init(&tsc);
    tsc_start(&tsc);
    
    int i; // , j, start, end;
    
    // start = (id * n) / ncpus;
    // end = ((id + 1) * n) / ncpus;
    
    for (i = 0; i < nloops; i++) {
        VECTOR_NAME(_init_from_map)(temp, 0, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
        spmv_mt_sym_fn(spm_mt_thread->spm, spm_mt_thread->data, y, y);
        pthread_barrier_wait(&barrier);
        /* for (j = 0; j < ncpus; j++)
            VECTOR_NAME(_add_part)(y, temp[j], y, start, end);
        */
	VECTOR_NAME(_add_from_map)(y, temp, y, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
    }
    
    tsc_pause(&tsc);
    secs = tsc_getsecs(&tsc);
    tsc_shut(&tsc);
    
    return NULL;
}

float SPMV_NAME(_bench_sym_mt_loop_numa)(spm_mt_t *spm_mt, unsigned long loops,
                                         unsigned long nrows,
                                         unsigned long ncols,
                                         SPMV_NAME(_fn_t) *fn)
{
	int err, i;
	pthread_t *tids;

	setaffinity_oncpu(spm_mt->spm_threads[0].cpu);
	
	nloops = loops;
	ncpus = spm_mt->nr_threads;
	assert(nrows == ncols);
	n = nrows;
	
	err = pthread_barrier_init(&barrier, NULL, ncpus);
	if (err){
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = (pthread_t *) malloc(ncpus * sizeof(pthread_t));
	if (!tids){
		perror("malloc");
		exit(1);
	}

	size_t *parts = malloc(ncpus * sizeof(*parts));
	int *nodes = (int *) malloc(ncpus * sizeof(int));
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	int nr_nodes = numa_max_node() + 1;
	VECTOR_TYPE **xs = malloc(sizeof(*xs)*nr_nodes);
	for (i = 0; i < nr_nodes; i++) {
		xs[i] = NULL;
	}

	VECTOR_TYPE *x_proto = xs[0];
	for (i = 0; i < ncpus; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		int node = spm_mt->spm_threads[i].node;
		if (!xs[node]) {
			xs[node] = VECTOR_NAME(_create_onnode)(n, node);
			if (!x_proto) {
				VECTOR_NAME(_init_rand_range)(xs[node],
				                              (ELEM_TYPE) -1000,
				                              (ELEM_TYPE) 1000);
				x_proto = xs[node];
			} else {
				/* copy the elements from the prototype */
				memcpy(xs[node]->elements, x_proto->elements,
				       n * sizeof(ELEM_TYPE));
			}
		}

		parts[i] = spm->nr_rows * sizeof(ELEM_TYPE);
		nodes[i] = spm->node;
		spm->data = xs[spm->node];
		if (fn)
			spm->spmv_fn = fn;
	}
	
	printf("check for allocation of x vector\n");
	for (i = 0; i < nr_nodes; i++)
		if (xs[i])
			check_onnode(xs[i]->elements,
			             n  * sizeof(*xs[i]->elements), i);

	/* Allocate an interleaved y */
	y = VECTOR_NAME(_create_interleaved)(n, parts, ncpus, nodes);
	VECTOR_NAME(_init)(y, 0);
	
	printf("check for allocation of y vector\n");
	check_interleaved(y->elements, n * sizeof(*y->elements), parts, ncpus,
	                  nodes);

	/* Allocate temporary buffers */
	temp = malloc(ncpus * sizeof(*temp));
	temp[0] = y;
	for (i = 1; i < ncpus; i++) {
		int tnode = spm_mt->spm_threads[i].node;
	    temp[i] = VECTOR_NAME(_create_onnode)(n, tnode);
	}
	
	for (i = 1; i < ncpus; i++) {
	    int j;
	    for (j = 1; j < n; j++)
	        temp[i]->elements[j] = 0;
	}
	
	printf("check for allocation of temp vectors\n");
	for (i = 1; i < ncpus; i++) {
		int tnode = spm_mt->spm_threads[i].node;
		
        	check_onnode(temp[i]->elements, n * sizeof(*temp[i]->elements),
		             tnode);
	}
    
	for (i = 1; i < ncpus; i++)
		pthread_create(tids + i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	do_spmv_thread_main(spm_mt->spm_threads);

	for (i = 1; i < ncpus; i++){
		pthread_join(tids[i], NULL);
	}

	/* Destroy vectors */
	for (i = 0; i < nr_nodes; i++) {
		if (xs[i])
			VECTOR_NAME(_destroy)(xs[i]);
	}
	free(xs);

	VECTOR_NAME(_destroy)(y);

    for (i = 1; i < ncpus; i++) {
		VECTOR_NAME(_destroy)(temp[i]);
	}
	free(temp);
	
	free(parts);
	free(nodes);
	free(tids);
	pthread_barrier_destroy(&barrier);
	return secs;
}

void SPMV_NAME(_check_sym_mt_loop_numa)(void *spm_serial, spm_mt_t *spm_mt,
                                        SPMV_NAME(_fn_t) *fn,
                                        unsigned long loops,
                                        unsigned long nrows,
                                        unsigned long ncols,
                                        SPMV_NAME(_fn_t) *mt_fn)
{
	int err, i;
	pthread_t *tids;
	VECTOR_TYPE *y2;

	nloops = loops;
	ncpus = spm_mt->nr_threads;
	assert(nrows == ncols);
	n = nrows;
	
	err = pthread_barrier_init(&barrier, NULL, ncpus + 1);
	if (err){
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = (pthread_t *) malloc(ncpus * sizeof(pthread_t));
	if ( !tids ){
		perror("malloc");
		exit(1);
	}

	size_t *parts = malloc(ncpus * sizeof(*parts));
	int *nodes = malloc(ncpus * sizeof(*nodes));
	if (!parts || !nodes) {
		perror("malloc");
		exit(1);
	}

	int nr_nodes = numa_max_node() + 1;
	VECTOR_TYPE **xs = malloc(sizeof(*xs)*nr_nodes);
	for (i = 0; i < nr_nodes; i++) {
		xs[i] = NULL;
	}

	VECTOR_TYPE *x_proto = NULL;
	for (i = 0; i < ncpus; i++) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		int node = spm->node;
		if (!xs[node]) {
			xs[node] = VECTOR_NAME(_create_onnode)(n, node);
			if (!x_proto) {
				VECTOR_NAME(_init_rand_range)(xs[node],
				                              (ELEM_TYPE) -1000,
				                              (ELEM_TYPE) 1000);
				x_proto = xs[node];
			} else {
				/* copy the elements from the prototype */
				memcpy(xs[node]->elements, x_proto->elements,
				       n * sizeof(ELEM_TYPE));
			}
		}

		parts[i] = spm->nr_rows * sizeof(ELEM_TYPE);
		nodes[i] = node;
		spm->data = xs[node];
		if (mt_fn)
			spm->spmv_fn = mt_fn;
	}
	
	/* Allocate an interleaved y */
	y = VECTOR_NAME(_create_interleaved)(n, parts, ncpus, nodes);
	y2 = VECTOR_NAME(_create)(n);
	VECTOR_NAME(_init)(y, 0);
	VECTOR_NAME(_init)(y2, 0);

	/* Allocate temporary buffers */
	temp = malloc(ncpus * sizeof(*temp));
	temp[0] = y;
	for (i = 1; i < ncpus; i++) {
		int node = spm_mt->spm_threads[i].node;
		
	    temp[i] = VECTOR_NAME(_create_onnode)(n, node);
	}
	
	for (i = 1; i < ncpus; i++) {
	    int j;
	    for (j = 1; j < n; j++)
	        temp[i]->elements[j] = 0;
	}
	
	for (i = 0; i < ncpus; i++)
		pthread_create(tids + i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (i = 0; i < nloops; i++) {
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		/* use the prototype of x */
		fn(spm_serial, x_proto, y2);
		if (VECTOR_NAME(_compare)(y2, y) < 0) {
			exit(1);
		}
	}

	/* Destroy vectors */
	for (i = 0; i < nr_nodes; i++) {
		if (xs[i])
			VECTOR_NAME(_destroy)(xs[i]);
	}
	free(xs);

	VECTOR_NAME(_destroy)(y);

    for (i = 1; i < ncpus; i++) {
		VECTOR_NAME(_destroy)(temp[i]);
	}
	free(temp);
	
	free(parts);
	free(nodes);
	free(tids);
	pthread_barrier_destroy(&barrier);
}
