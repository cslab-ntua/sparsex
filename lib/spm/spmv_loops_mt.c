/*
 * spmv_loops_mt.c
 *
 * Copyright (C) 2007-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <stdlib.h>
#include <pthread.h>

#include "spmv_method.h"
#include "vector.h"
#include "spm_mt.h"
#include "mt_lib.h"

#include "spmv_loops_mt.h"
#include "tsc.h"
#ifdef SPMV_PRFCNT
#include "prfcnt.h"
#endif /* SPMV_PRFCNT */

static VECTOR_TYPE *x = NULL;
static VECTOR_TYPE *y = NULL;
static pthread_barrier_t barrier;
static unsigned long loops_nr = 0;
static float secs = 0.0;

static void *do_spmv_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
	SPMV_NAME(_fn_t) *spmv_mt_fn = spm_mt_thread->spmv_fn;
	setaffinity_oncpu(spm_mt_thread->cpu);
	int i;
    tsc_t thread_tsc;

#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt = (prfcnt_t *) spm_mt_thread->data;
	prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif

    tsc_init(&thread_tsc);
	for (i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
        tsc_start(&thread_tsc);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
        tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
	}

#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif

    spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
    tsc_shut(&thread_tsc);
	return NULL;
}

#define SWAP(x,y) \
	do { \
		typeof(x) _tmp; \
		_tmp = x; \
		x = y; \
		y = _tmp; \
	} while (0)

static void *do_spmv_thread_main_swap(void *arg)
{
	spm_mt_thread_t *spm_mt_thread;
#ifdef SPMV_PRFCNT
	prfcnt_t *prfcnt;
#endif
	SPMV_NAME(_fn_t) *spmv_mt_fn;
	tsc_t total_tsc, thread_tsc;

	spm_mt_thread = arg;
	spmv_mt_fn = spm_mt_thread->spmv_fn;
#ifdef SPMV_PRFCNT
	prfcnt = (prfcnt_t *) spm_mt_thread->data;
#endif
	setaffinity_oncpu(spm_mt_thread->cpu);

	VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);

	// Assert this is a square matrix and swap is ok.
	assert(x->size == y->size);
	tsc_init(&total_tsc);
	tsc_init(&thread_tsc);
	tsc_start(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
	prfcnt_start(prfcnt);
#endif
	int i;
	for (i = 0; i < loops_nr; i++) {
		pthread_barrier_wait(&barrier);
		tsc_start(&thread_tsc);
		spmv_mt_fn(spm_mt_thread->spm, x, y);
		tsc_pause(&thread_tsc);
		pthread_barrier_wait(&barrier);
		SWAP(x, y);
	}
	tsc_pause(&total_tsc);
#ifdef SPMV_PRFCNT
	prfcnt_pause(prfcnt);
#endif
	spm_mt_thread->secs = tsc_getsecs(&thread_tsc);
	secs = tsc_getsecs(&total_tsc);
	tsc_shut(&total_tsc);
	tsc_shut(&thread_tsc);

	return NULL;
}

float SPMV_NAME(_bench_mt_loop)(spm_mt_t *spm_mt, unsigned long loops,
                                unsigned long rows_nr, unsigned long cols_nr,
                                SPMV_NAME(_fn_t) *fn)
{
	int i, err;
	pthread_t *tids;

	secs = 0.0;
	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);
	loops_nr = loops;

	err = pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads);
	if (err) {
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = malloc(sizeof(pthread_t)*spm_mt->nr_threads);
	if (!tids) {
		perror("malloc");
		exit(1);
	}

	for (i = 0; i < spm_mt->nr_threads; i++) {
		if (fn != NULL)
			spm_mt->spm_threads[i].spmv_fn = fn;
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
	}

	/*
	 * spawn two kind of threads:
	 *	- 1	: do_spmv_thread_main_swap -> computes and does swap.
	 *	- N-1	: do_spmv_thread -> just computes.
	 */
	pthread_create(tids, NULL, do_spmv_thread_main_swap, spm_mt->spm_threads);
	for (i = 1; i < spm_mt->nr_threads; i++)
		pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);

	for (i = 0; i < spm_mt->nr_threads; i++)
		pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
	// Report performance counters for every thread.
	for (i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		fprintf(stdout, "Perf. counters: thread %d on cpu %d\n", i,
		        spmv_thread.cpu);
		prfcnt_report(prfcnt);
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	pthread_barrier_destroy(&barrier);
	free(tids);
	return secs;
}

void SPMV_NAME(_check_mt_loop)(void *spm, spm_mt_t *spm_mt,
                               SPMV_NAME(_fn_t) *fn, unsigned long loops,
                               unsigned long rows_nr, unsigned long cols_nr,
                               SPMV_NAME(_fn_t) *mt_fn)
{
	int err, i;
	pthread_t *tids;
	VECTOR_TYPE *y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);
	y2 = VECTOR_NAME(_create)(rows_nr);
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

	for (i = 0; i < spm_mt->nr_threads; i++) {
		if (mt_fn != NULL)
			spm_mt->spm_threads[i].spmv_fn = mt_fn;
#ifdef SPMV_PRFCNT
		spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
		if (!spm_mt->spm_threads[i].data) {
			perror("malloc");
			exit(1);
		}
#endif
		pthread_create(tids+i, NULL, do_spmv_thread,
		               spm_mt->spm_threads + i);
	}

	for (i = 0; i < loops; i++) {
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);
		VECTOR_NAME(_init)(y2, (ELEM_TYPE) 21);
		pthread_barrier_wait(&barrier);
		pthread_barrier_wait(&barrier);
		fn(spm, x, y2);
		if (VECTOR_NAME(_compare)(y2, y) < 0)
			exit(1);
	}

	for (i = 0; i < spm_mt->nr_threads; i++)
		pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
	for (i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
		prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
		prfcnt_shut(prfcnt);
		free(prfcnt);
	}
#endif

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	VECTOR_NAME(_destroy)(y2);
	pthread_barrier_destroy(&barrier);
	free(tids);
}

void SPMV_NAME(_check_mt_loop_serial)(void *spm, spm_mt_t *spm_mt,
                                      SPMV_NAME(_fn_t) *fn,
                                      unsigned long loops,
                                      unsigned long rows_nr,
                                      unsigned long cols_nr,
                                      SPMV_NAME(_fn_t) *mt_fn)
{
	int i, j;
	VECTOR_TYPE *y2;

	x = VECTOR_NAME(_create)(cols_nr);
	y = VECTOR_NAME(_create)(rows_nr);
	y2 = VECTOR_NAME(_create)(rows_nr);
	loops_nr = loops;

	for (i = 0; i < spm_mt->nr_threads; i++)
		if (mt_fn != NULL)
			spm_mt->spm_threads[i].spmv_fn = mt_fn;

	for (i = 0; i < loops; i++) {
		VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);
		VECTOR_NAME(_init)(y, (ELEM_TYPE) 0);
		VECTOR_NAME(_init)(y2, (ELEM_TYPE) 21);

		for (j = 0; j < spm_mt->nr_threads; j++) {
			spm_mt_thread_t *t = spm_mt->spm_threads + j;
			SPMV_NAME(_fn_t) *spmv_mt_fn = t->spmv_fn;
			spmv_mt_fn(t->spm, x, y);
		}

		fn(spm, x, y2);
		if (VECTOR_NAME(_compare)(y2, y) < 0)
			exit(1);
	}

	VECTOR_NAME(_destroy)(x);
	VECTOR_NAME(_destroy)(y);
	VECTOR_NAME(_destroy)(y2);
}
