/*
 * spmv_matvec_mt.c -- Multithread matvec product routines.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011, Vasileios Karakasis
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
#include "spmv_matvec_mt.h"

static pthread_barrier_t barrier;

static void *do_matvec_thread(void *arg)
{
	spm_mt_thread_t *spm_mt_thread = (spm_mt_thread_t *) arg;
	SPMV_NAME(_fn_t) *spmv_mt_fn = spm_mt_thread->spmv_fn;

	setaffinity_oncpu(spm_mt_thread->cpu);
	pthread_barrier_wait(&barrier);
#ifndef _CSR_
    VECTOR_NAME(_init_part)(spm_mt_thread->y, spm_mt_thread->row_start,
                            spm_mt_thread->nr_rows, (ELEM_TYPE) 0);
#endif
	spmv_mt_fn(spm_mt_thread->spm, spm_mt_thread->x, spm_mt_thread->y);
	pthread_barrier_wait(&barrier);
	return (void *) 0;
}

void SPMV_NAME(_matvec_mt)(spm_mt_t *spm_mt, VECTOR_TYPE *x, VECTOR_TYPE *y)
{
	pthread_t *tids;
	int i;

	if (pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads)) {
		perror("pthread_barrier_init");
		exit(1);
	}

	tids = malloc(sizeof(*tids)*spm_mt->nr_threads);
    if (!tids) {
        perror("malloc() failed");
        exit(1);
    }

	for (i = 0; i < spm_mt->nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
	}

	for (i = 1; i < spm_mt->nr_threads; i++) {
		if (pthread_create(tids + i, NULL, do_matvec_thread,
						   spm_mt->spm_threads + i)) {
			perror("pthread_create");
			exit(1);
		}
	}

    do_matvec_thread(spm_mt->spm_threads);

	for (i = 1; i < spm_mt->nr_threads; i++) {
		if (pthread_join(tids[i], NULL)) {
			perror("pthread_join");
			exit(1);
		}
	}

	free(tids);
}
