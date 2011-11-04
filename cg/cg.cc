/*
 * cg.cc -- The CG Manager Implementation.
 *
 * Copyright (C) 2011,      Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
 
#include "cg.h"

void FindSolution(spm_mt_t *spm_mt, vector_double_t *sol, vector_double_t *b, 
                  vector_double_t *temp)
{
    unsigned long i;
    spmv_double_fn_t *fn;
    spm_mt_thread_t *spm_thread;
    
    vector_double_init(b, 0);
    for (i = 0; i < spm_mt->nr_threads; i++) {
        vector_double_init(temp, 0);
        spm_thread = spm_mt->spm_threads + i;
        fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
        fn(spm_thread->spm, sol, temp);
        vector_double_add(b, temp, b);
    }
}

void FindSymSolution(spm_mt_t *spm_mt, vector_double_t *sol,
                     vector_double_t *b, vector_double_t *temp)
{
    unsigned long i;
    spmv_double_sym_fn_t *fn;
    spm_mt_thread_t *spm_thread;
    
    vector_double_init(b, 0);
    for (i = 0; i < spm_mt->nr_threads; i++) {
        vector_double_init(temp, 0);
        spm_thread = spm_mt->spm_threads + i;
        fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
        fn(spm_thread->spm, sol, temp, temp);
        vector_double_add(b, temp, b);
    }
}

void InitializeCg(spm_mt_t *spm_mt, vector_double_t *x, vector_double_t *r,
                  vector_double_t *p, vector_double_t *b, vector_double_t *temp)
{
    unsigned long i;
    spmv_double_fn_t *fn;
    spm_mt_thread_t *spm_thread;

    vector_double_init(x, 0.01);                        ///> Initialize x
                                                        //   to 0.01.
    vector_double_init(temp, 0);                        ///> Do temp = A * x0.
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
        vector_double_init(r, 0);
        fn(spm_thread->spm, x, r);
        vector_double_add(temp, r, temp);
    }
    
    vector_double_sub(b, temp, r);                      ///> Do r0 = b - temp.
    vector_double_copy(r, p);                           ///> Do p0 = r0.

}

void InitializeSymCg(spm_mt_t *spm_mt, vector_double_t *x, vector_double_t *r,
                     vector_double_t *p, vector_double_t *b,
                     vector_double_t *temp)
{
    unsigned long i;
    spmv_double_sym_fn_t *fn;
    spm_mt_thread_t *spm_thread;

    vector_double_init(x, 0.01);                        ///> Initialize x
                                                        //   to 0.01.
    vector_double_init(temp, 0);                        ///> Do temp = A * x0.
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
        vector_double_init(r, 0);
        fn(spm_thread->spm, x, r, r);
        vector_double_add(temp, r, temp);
    }
    
    vector_double_sub(b, temp, r);                      ///> Do r0 = b - temp.
    vector_double_copy(r, p);                           ///> Do p0 = r0.

}


void *NormalCgSideThread(void *arg)
{
    uint64_t i;
    cg_params *params = (cg_params *) arg;
    uint64_t nloops = params->nloops;
    spm_mt_thread_t *spm_thread = params->spm_thread;
    spmv_double_fn_t *fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
    vector_double_t *p = params->p;
    vector_double_t *temp = params->temp;
    vector_double_t *r = params->r;
    vector_double_t *x = params->x;
    double *rr = params->rr;
    double *tp = params->tp;
    double *rr_new = params->rr_new;
    double *ai = params->ai;
    double *bi = params->bi;
    uint64_t start = params->start;
    uint64_t end = params->end;
    pthread_barrier_t *barrier = params->barrier;
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(spm_thread->cpu);
    
    ///> Do nr_loops SpMxVs.
    for (i = 0; i < nloops; i++) {
        pthread_barrier_wait(barrier);
        fn(spm_thread->spm, p, temp);
        pthread_barrier_wait(barrier);
        *rr = vector_double_mul_part(r, r, start, end);
        *tp = vector_double_mul_part(temp, p, start, end);
        pthread_barrier_wait(barrier);
        pthread_barrier_wait(barrier);
        vector_double_scale_add_part(r, temp, r, -(*ai), start, end);
        vector_double_scale_add_part(x, p, x, *ai, start, end);
        pthread_barrier_wait(barrier);
        *rr_new = vector_double_mul_part(r, r, start, end);
        pthread_barrier_wait(barrier);
        pthread_barrier_wait(barrier);
	vector_double_scale_add_part(r, p, p, *bi, start, end);
        pthread_barrier_wait(barrier);
    }
    
    return NULL;
}

void *SymCgSideThread(void *arg)
{
    uint64_t i;
    cg_params *params = (cg_params *) arg;
    uint64_t nloops = params->nloops;
    spm_mt_thread_t *spm_thread = params->spm_thread;
    uint64_t id = spm_thread->id;
    spmv_double_sym_fn_t *fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
    vector_double_t *p = params->p;
    vector_double_t *temp = params->temp;
    vector_double_t *r = params->r;
    vector_double_t *x = params->x;
    double *rr = params->rr;
    double *tp = params->tp;
    double *rr_new = params->rr_new;
    double *ai = params->ai;
    double *bi = params->bi;
    pthread_barrier_t *barrier = params->barrier;
    vector_double_t **sub_vectors = params->sub_vectors;
    uint64_t start = params->start;
    uint64_t end = params->end;
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(spm_thread->cpu);
    
    ///> Do nr_loops SpMxVs.
    for (i = 0; i < nloops; i++) {
        pthread_barrier_wait(barrier);
        fn(spm_thread->spm, p, temp, sub_vectors[id]);
        pthread_barrier_wait(barrier);
        vector_double_add_from_map(temp, sub_vectors, temp, spm_thread->map);
        vector_double_init_from_map(sub_vectors, 0, spm_thread->map);
        pthread_barrier_wait(barrier);
        *rr = vector_double_mul_part(r, r, start, end);
        *tp = vector_double_mul_part(temp, p, start, end);
        pthread_barrier_wait(barrier);
        pthread_barrier_wait(barrier);
        vector_double_scale_add_part(r, temp, r, -(*ai), start, end);
        vector_double_scale_add_part(x, p, x, *ai, start, end);
        pthread_barrier_wait(barrier);
        *rr_new = vector_double_mul_part(r, r, start, end);
        pthread_barrier_wait(barrier);
        pthread_barrier_wait(barrier);
	vector_double_scale_add_part(r, p, p, *bi, start, end);
        pthread_barrier_wait(barrier);
    }
    
    return NULL;
}

void NormalCgMainThread(cg_params *params, double *cg_time, double *spmv_time,
                        double * red_time)
{
    uint64_t i, j;
    uint64_t nloops = params->nloops;
    uint64_t ncpus = params->ncpus;
    spm_mt_thread_t *spm_thread = params->spm_thread;
    spmv_double_fn_t *fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
    vector_double_t *p = params->p;
    vector_double_t *temp = params->temp;
    vector_double_t *r = params->r;
    vector_double_t *x = params->x;
    double *rr = params->rr;
    double *tp = params->tp;
    double *rr_new = params->rr_new;
    double *ai = params->ai;
    double *bi = params->bi;
    uint64_t start = params->start;
    uint64_t end = params->end;
    pthread_barrier_t *barrier = params->barrier;
    
    // Set thread to the appropriate cpu.
    setaffinity_oncpu(spm_thread->cpu);
    
    // Initialize timers.
    timer_init(&cg_timer);
    timer_init(&spmv_timer);
    timer_init(&reduction_timer);
    
    // Start main timer.
    timer_start(&cg_timer);
    
    ///> Do nr_loops.
    for (i = 0; i < nloops; i++) {
        pthread_barrier_wait(barrier);
        timer_start(&spmv_timer);
        
        ///> Do temp = A*p.
        fn(spm_thread->spm, p, temp);
        pthread_barrier_wait(barrier);
        timer_pause(&spmv_timer);
        
        ///> Calculate ai.
        rr[0] = vector_double_mul_part(r, r, start, end);
        tp[0] = vector_double_mul_part(temp, p, start, end);
        pthread_barrier_wait(barrier);
        for (j = 1; j < ncpus; j++) {
            rr[0] += rr[j];
            tp[0] += tp[j];
        }
        *ai = rr[0] / tp[0];
        
        pthread_barrier_wait(barrier);
        ///> Do r = r - ai*A*p.
        vector_double_scale_add_part(r, temp, r, -(*ai), start, end);
        ///> Do x = x + ai*p.
        vector_double_scale_add_part(x, p, x, *ai, start, end);
        pthread_barrier_wait(barrier);
        
        ///> Calculate bi.
	rr_new[0] = vector_double_mul_part(r, r, start, end);
        pthread_barrier_wait(barrier);
        for (j = 1; j < ncpus; j++)
            rr_new[0] += rr_new[j];
        *bi = rr_new[0] / rr[0];
        pthread_barrier_wait(barrier);
        
        ///> Do p = r + bi*p.
	vector_double_scale_add_part(r, p, p, *bi, start, end);
        pthread_barrier_wait(barrier);
    }
    
    timer_pause(&cg_timer);
    
    *cg_time = timer_secs(&cg_timer);
    *spmv_time = timer_secs(&spmv_timer);
    *red_time = timer_secs(&reduction_timer);
}

void SymCgMainThread(cg_params *params, double *cg_time, double *spmv_time,
                     double * red_time)
{
    uint64_t i, j;
    uint64_t nloops = params->nloops;
    uint64_t ncpus = params->ncpus;
    spm_mt_thread_t *spm_thread = params->spm_thread;
    uint64_t id = spm_thread->id;
    spmv_double_sym_fn_t *fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
    vector_double_t *p = params->p;
    vector_double_t *temp = params->temp;
    vector_double_t *r = params->r;
    vector_double_t *x = params->x;
    double *rr = params->rr;
    double *tp = params->tp;
    double *rr_new = params->rr_new;
    double *ai = params->ai;
    double *bi = params->bi;
    pthread_barrier_t *barrier = params->barrier;
    vector_double_t **sub_vectors = params->sub_vectors;  
    uint64_t start = params->start;
    uint64_t end = params->end;
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(spm_thread->cpu);
    
    // Initialize timers.
    timer_init(&cg_timer);
    timer_init(&spmv_timer);
    timer_init(&reduction_timer);
    
    // Start main timer.
    timer_start(&cg_timer);
    
    ///> Do nr_loops.
    for (i = 0; i < nloops; i++) {
        pthread_barrier_wait(barrier);
        timer_start(&spmv_timer);
        fn(spm_thread->spm, p, temp, sub_vectors[id]);  ///> Do temp = Ap.
        pthread_barrier_wait(barrier);
        timer_pause(&spmv_timer);
        timer_start(&reduction_timer);
        vector_double_add_from_map(temp, sub_vectors, temp, spm_thread->map);
        vector_double_init_from_map(sub_vectors, 0, spm_thread->map);
        pthread_barrier_wait(barrier);
        timer_pause(&reduction_timer);
        
        ///> Calculate ai.
        rr[0] = vector_double_mul_part(r, r, start, end);
        tp[0] = vector_double_mul_part(temp, p, start, end);
        pthread_barrier_wait(barrier);
        for (j = 1; j < ncpus; j++) {
            rr[0] += rr[j];
            tp[0] += tp[j];
        }
        *ai = rr[0] / tp[0];
        
        pthread_barrier_wait(barrier);
        ///> Do r = r - ai*A*p.
        vector_double_scale_add_part(r, temp, r, -(*ai), start, end);
        ///> Do x = x + ai*p.
        vector_double_scale_add_part(x, p, x, *ai, start, end);
        pthread_barrier_wait(barrier);
        
        ///> Calculate bi.
	rr_new[0] = vector_double_mul_part(r, r, start, end);
        pthread_barrier_wait(barrier);
        for (j = 1; j < ncpus; j++)
            rr_new[0] += rr_new[j];
        *bi = rr_new[0] / rr[0];
        pthread_barrier_wait(barrier);
        
        ///> Do p = r + bi*p.
	vector_double_scale_add_part(r, p, p, *bi, start, end);
        pthread_barrier_wait(barrier);
    }
    
    timer_pause(&cg_timer);
    
    *cg_time = timer_secs(&cg_timer);
    *spmv_time = timer_secs(&spmv_timer);
    *red_time = timer_secs(&reduction_timer);
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
