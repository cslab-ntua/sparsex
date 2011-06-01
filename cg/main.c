/* -*- C -*-
 *
 * main.c -- Implementation of Conjugate Gradient Method (CG).
 *
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cblas.h>
#include <math.h>
#include <libgen.h>
#include <assert.h>
#include <pthread.h>
#include <sys/time.h>

#include "cg_vector.h"
#include "cg.h"

///> Type of method used in CG.
typedef enum {
   CSR = 0,
   CSR_s = 1,
   CSX = 2,
   CSX_s = 3
} cg_method_t;


int main(int argc, char *argv[])
{
    unsigned long       i;
    char                *mmf_file;
    uint64_t            nr_rows, nr_cols, nr_nzeros;
    spm_mt_t            *spm_mt;
    spm_mt_thread_t     *spm_thread;
    vector_double_t     *x, *b, *sol, *r, *p, *t;
    double              rr, tp, ai, bi, rr_new;
    struct timeval      start, end;
    double              time, avg_error, sol_dis;
    pthread_t           *tids;
    spmv_double_fn_t    *fn;
    
    ///> Default options.
    unsigned long   nr_loops = 128;
    cg_method_t     cg_method = CSR;
    
    ///> Parse Options.
    while ((i = getopt(argc, argv, "xsl:")) != -1) {
        switch (i) {
            case 'x':
                cg_method = CSX;
                break;
            
            case 's':
                cg_method++;
                break;
            
            case 'l':
                nr_loops = atol(optarg);
                break;
                
	    default:
	        fprintf(stderr, "Usage: cg [-x -l <number of loops>] mmf_file\n");
	        exit(1);
	}
    }
    
    ///> Take input file from which matrix will be loaded.
    if (argc != optind+1) {
        fprintf(stderr, "Usage: cg [-x -l <number of loops>] mmf_file\n");
        exit(1);
    }
    mmf_file = argv[optind];
    
    ///> Load matrix in appropriate format.
    switch(cg_method) {
        case(CSR) :
            spm_mt = (spm_mt_t *) spm_crs32_double_mt_init_mmf(mmf_file, &nr_rows, &nr_cols, &nr_nzeros);
	    method_t *meth = method_get("spm_crs32_double_mt_multiply");
            for (i = 0; i < spm_mt->nr_threads; i++)
                spm_mt->spm_threads[i].spmv_fn = meth->fn;
            break;
        default :
            fprintf(stderr, "Wrong method choosed\n");
            exit(1);
    }
            
    
    ///> Init pthreads.
    tids = (pthread_t *) malloc((spm_mt->nr_threads - 1) * sizeof(pthread_t));
    if (!tids) {
	fprintf(stderr, "Malloc of pthreads failed\n");
	exit(1);
    }
    if (pthread_barrier_init(&barrier, NULL, spm_mt->nr_threads)) {
        fprintf(stderr, "Pthread barrier init failed");
        exit(1);
    }

    ///> Create vectors with the appropriate size.
    x = vector_double_create(nr_cols);
    sol = vector_double_create(nr_cols);
    b = vector_double_create(nr_rows);
    r = vector_double_create(nr_cols);
    p = vector_double_create(nr_cols);
    t = vector_double_create(nr_cols);
    
    ///> Assign the appropriate parameters to each thread.
    cg_params *params = (cg_params *) malloc(spm_mt->nr_threads * sizeof(cg_params));
    for (i = 0; i < spm_mt->nr_threads; i++) {
        params[i].nr_loops = nr_loops;
        params[i].spm_thread = spm_mt->spm_threads+i;
        params[i].in = p;
        params[i].out = t;
    }
    
    ///> Load vector sol with random values.
    vector_double_init_rand_range(sol, (double) -1000, (double) 1000);
	
    ///> Do A*sol=b to define b vector.
    vector_double_init(b, 0);
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        fn = spm_thread->spmv_fn;
        vector_double_init(t, 0);
        fn(spm_thread->spm, sol, t);
        vector_double_add(b, t, b);
    }
    sol_dis = cblas_ddot(nr_cols, sol->elements, 1, sol->elements, 1);
    sol_dis = sqrt(sol_dis);
    
    vector_double_init(x, 0.01);                                ///> Initiate vector x (x0).
    vector_double_init(t, 0);                                   ///> Do t = A * x0.
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        fn = spm_thread->spmv_fn;
        vector_double_init(r, 0);
        fn(spm_thread->spm, x, r);
        vector_double_add(t, r, t);
    }
    vector_double_sub(b, t, r);                                 ///> Do r0 = b - t.
    vector_double_copy(r, p);                                   ///> Do p0 = r0.

    ///> Start counting time.
    gettimeofday(&start, NULL);
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(params[0].spm_thread->cpu);
    
    ///> Initiate side threads.
    for (i = 1; i < spm_mt->nr_threads; i++)
	pthread_create(&tids[i-1], NULL, cg_side_thread, (void *) &params[i]);

    ///> Execute main thread.
    spm_thread = params[0].spm_thread;
    fn = spm_thread->spmv_fn;
    for (i = 0; i < nr_loops; i++) {
        pthread_barrier_wait(&barrier);
        fn(spm_thread->spm, p, t);                                      ///> Do t = Ap.
        pthread_barrier_wait(&barrier);
        rr = cblas_ddot(nr_cols, r->elements, 1, r->elements, 1);       ///> Calculate ai.
        tp = cblas_ddot(nr_cols, t->elements, 1, p->elements, 1);
        ai = rr/tp;
        cblas_daxpy(nr_cols, ai, p->elements, 1, x->elements, 1);       ///> Do x = x + ai*p.
	cblas_daxpy(nr_cols, -ai, t->elements, 1, r->elements, 1);      ///> Do r = r - ai*A*p.
	/*if (i % 10 == 0) {
            vector_double_sub(sol, x, t);
            avg_error = cblas_ddot(nr_cols, t->elements, 1, t->elements, 1);
            avg_error = sqrt(avg_error);
            avg_error /= sol_dis;
            printf("Loop: %lu Relative Distance: %lf\n", i, avg_error);
	}*/
	rr_new = cblas_ddot(nr_cols, r->elements, 1, r->elements, 1);   ///> Calculate bi.
	bi = rr_new / rr;
	vector_double_init(t, 0);                                       ///> Do p = r + bi*p.
	cblas_daxpy(nr_cols, bi, p->elements, 1, t->elements, 1);
        vector_double_add(r, t, p);
    }
    
    ///> Stop counting time.
    gettimeofday(&end, NULL);
    time = end.tv_sec-start.tv_sec + (end.tv_usec-start.tv_usec)/1000000.0;
    
    ///> Print results.
    vector_double_sub(sol, x, t);
    avg_error = cblas_ddot(nr_cols, t->elements, 1, t->elements, 1);
    avg_error = sqrt(avg_error);
    avg_error /= sol_dis;
    //printf("Loop: %lu Relative Distance: %lf\n", nr_loops, avg_error);
    printf("m:%s l:%lu rd:%lf t:%lf\n", basename(mmf_file), nr_loops, avg_error, time);
    
    ///> Release vectors.
    vector_double_destroy(x);
    vector_double_destroy(sol);
    vector_double_destroy(b);
    vector_double_destroy(r);
    vector_double_destroy(p);
    vector_double_destroy(t);
    
    ///> Free pthreads.
    free(tids);
    
    ///> Release matrix.
    // TODO:Destroy Matrix.
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
