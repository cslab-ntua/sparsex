/*
 * spmv_loops_sym_mt.c
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "spmv_loops_sym_mt.h"

static VECTOR_TYPE *x = NULL;
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
    int i, start, end;
    
    setaffinity_oncpu(spm_mt_thread->cpu);

#ifdef SPMV_PRFCNT
    prfcnt_t *prfcnt = (prfcnt_t *) spm_mt_thread->data;
    prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
    prfcnt_start(prfcnt);
#endif

    start = (id * n) / ncpus;
    end = ((id + 1) * n) / ncpus;
    
    for (i = 0; i < nloops; i++) {
        /*if (id != 0)
            VECTOR_NAME(_init)(temp[id], 0);*/
        VECTOR_NAME(_init_from_map)(temp, 0, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
        spmv_mt_sym_fn(spm_mt_thread->spm, x, y, temp[id]);
        pthread_barrier_wait(&barrier);
        /*for (j = 1; j < ncpus; j++)
            VECTOR_NAME(_add_part)(y, temp[j], y, start, end);*/
        VECTOR_NAME(_add_from_map)(y, temp, y, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
    }
    
#ifdef SPMV_PRFCNT
    prfcnt_pause(prfcnt);
#endif

    return NULL;
}

#define SWAP(x,y)        \
    do {                 \
        typeof(x) _tmp;  \
        _tmp = x;        \
        x = y;           \
        y = _tmp;        \
    } while (0)

static void *do_spmv_thread_main_swap(void *arg)
{
    spm_mt_thread_t *spm_mt_thread;

#ifdef SPMV_PRFCNT
    prfcnt_t *prfcnt;
#endif

    SPMV_NAME(_sym_fn_t) *spmv_mt_sym_fn;
    int id;
    tsc_t tsc;
    
    spm_mt_thread = arg;
    spmv_mt_sym_fn = spm_mt_thread->spmv_fn;
    
    id = spm_mt_thread->id;

#ifdef SPMV_PRFCNT
    prfcnt = (prfcnt_t *) spm_mt_thread->data;
#endif
    setaffinity_oncpu(spm_mt_thread->cpu);

    VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE)-1000, (ELEM_TYPE)1000);

    // assert that this is a rectangular matrix, and swap is OK
    assert(x->size == y->size);
    tsc_init(&tsc);
    tsc_start(&tsc);
    
#ifdef SPMV_PRFCNT
    prfcnt_init(prfcnt, spm_mt_thread->cpu, PRFCNT_FL_T0 | PRFCNT_FL_T1);
    prfcnt_start(prfcnt);
#endif

    int i, start, end;
    
    start = (id * n) / ncpus;
    end = ((id + 1) * n) / ncpus;
    
    for (i = 0; i < nloops; i++) {
        VECTOR_NAME(_init_from_map)(temp, 0, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
        spmv_mt_sym_fn(spm_mt_thread->spm, x, y, y);
        pthread_barrier_wait(&barrier);
        /*for (j = 0; j < ncpus; j++)
            VECTOR_NAME(_add_part)(y, temp[j], y, start, end);*/
        VECTOR_NAME(_add_from_map)(y, temp, y, spm_mt_thread->map);
        pthread_barrier_wait(&barrier);
        SWAP(x, y);
    }
    tsc_pause(&tsc);

#ifdef SPMV_PRFCNT
    prfcnt_pause(prfcnt);
#endif

    secs = tsc_getsecs(&tsc);
    tsc_shut(&tsc);
    
    return NULL;
}

static void init_barrier(unsigned count)
{
    int err;
    
    err = pthread_barrier_init(&barrier, NULL, count);
    if (err) {
        perror("pthread_barrier_init");
        exit(1);
    }
}

float SPMV_NAME(_bench_sym_mt_loop)(spm_mt_t *spm_mt, unsigned long loops,
                                    unsigned long nrows, unsigned long ncols,
                                    SPMV_NAME(_fn_t) *fn)
{
    int i;
    pthread_t *tids;

    secs = 0.0;
    nloops = loops;
    ncpus = spm_mt->nr_threads;
    n = nrows;
    
    x = VECTOR_NAME(_create)(n);
    y = VECTOR_NAME(_create)(n);
    temp = (VECTOR_TYPE **) malloc(ncpus * sizeof(VECTOR_TYPE *));
    temp[0] = y;
    for (i = 1; i < ncpus; i++)
        temp[i] = VECTOR_NAME(_create)(n);

    init_barrier(ncpus);
    tids = (pthread_t *) malloc(ncpus * sizeof(pthread_t));
    if (!tids) {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }
    
    for (i = 0; i < ncpus; i++) {
        if (fn != NULL)
            spm_mt->spm_threads[i].spmv_fn = fn;

#ifdef SPMV_PRFCNT
        spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
        if (!spm_mt->spm_threads[i].data) {
            fprintf(stderr, "malloc failed\n");
            exit(1);
        }
#endif

    }
    
    pthread_create(tids, NULL, do_spmv_thread_main_swap, spm_mt->spm_threads);
    for (i = 1; i < ncpus; i++)
        pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
        
    for (i = 0; i < ncpus; i++)
        pthread_join(tids[i], NULL);

#ifdef SPMV_PRFCNT
    /* Report performance counters for every thread */
    for (i = 0; i < ncpus; i++) {
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
    for (i = 1; i < ncpus; i++)
        VECTOR_NAME(_destroy)(temp[i]);
    free(temp);
    free(tids);

    return secs;
}

void SPMV_NAME(_check_sym_mt_loop) (void *spm, spm_mt_t *spm_mt,
                                    SPMV_NAME(_fn_t) *fn, unsigned long loops,
                                    unsigned long nrows, unsigned long ncols,
                                    SPMV_NAME(_fn_t) *mt_fn)
{
    int i, j;
    pthread_t *tids;
    VECTOR_TYPE *y2;
    
    nloops = loops;
    ncpus = spm_mt->nr_threads;
    n = nrows;
    
    x = VECTOR_NAME(_create)(n);
    y = VECTOR_NAME(_create)(n);
    y2 = VECTOR_NAME(_create)(n);
    temp = (VECTOR_TYPE **) malloc(ncpus * sizeof(VECTOR_TYPE *));
    for (i = 0; i < ncpus; i++)
        temp[i] = VECTOR_NAME(_create)(n);

    init_barrier(ncpus + 1);
    tids = (pthread_t *) malloc(ncpus * sizeof(pthread_t));
    if (!tids) {
        fprintf(stderr, "malloc failed");
        exit(1);
    }

    for (i = 0; i < ncpus; i++) {
        if (mt_fn != NULL)
            spm_mt->spm_threads[i].spmv_fn = mt_fn;
            
#ifdef SPMV_PRFCNT
        spm_mt->spm_threads[i].data = malloc(sizeof(prfcnt_t));
        if (!spm_mt->spm_threads[i].data) {
            fprintf(stderr, "malloc failed");
            exit(1);
        }
#endif
        pthread_create(tids+i, NULL, do_spmv_thread, spm_mt->spm_threads + i);
    }

    for (i = 0; i < nloops; i++) {
        VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -1000, (ELEM_TYPE) 1000);
        VECTOR_NAME(_init)(y, (ELEM_TYPE) 21);
        VECTOR_NAME(_init)(y2, (ELEM_TYPE) 22);
        for (j = 1; j < ncpus; j++)
            VECTOR_NAME(_init)(temp[j], (ELEM_TYPE) 0);
        pthread_barrier_wait(&barrier);
        pthread_barrier_wait(&barrier);
        pthread_barrier_wait(&barrier);
        fn(spm, x, y2);
        if (VECTOR_NAME(_compare)(y2, y) < 0)
            exit(1);
    }

    for (i = 0; i < spm_mt->nr_threads; i++)
        pthread_join(tids[i], NULL);

    VECTOR_NAME(_destroy)(x);
    VECTOR_NAME(_destroy)(y);
    VECTOR_NAME(_destroy)(y2);
    for (i = 0; i < ncpus; i++)
        VECTOR_NAME(_destroy)(temp[i]);

#ifdef SPMV_PRFCNT
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_mt_thread_t spmv_thread = spm_mt->spm_threads[i];
        prfcnt_t *prfcnt = (prfcnt_t *) spmv_thread.data;
        prfcnt_shut(prfcnt);
        free(prfcnt);
    }
#endif

    free(tids);
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
