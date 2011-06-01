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

int main(int argc, char *argv[])
{
    unsigned long   i;
    char            *mmf_file;
    uint64_t        nr_rows, nr_cols, nr_nzeros, nr_loops;
    spm_mt_t        *A;
    vector_double_t *x, *b, *sol, *r, *p, *t;
    double          rr, tp, ai, bi, rr_new;
    struct timeval  start, end;
    double          time, avg_error, sol_dis;
    pthread_t       *tids;
    int             err;
    
    ///> Take input file from which matrix will be loaded.
    assert(argc == 3 && "usage: ./cg <number of loops> <filename of matrix>");
    mmf_file = argv[2];
    
    ///> Load matrix in CSR format.
    A = (spm_mt_t *) spm_crs32_double_mt_init_mmf(mmf_file, &nr_rows, &nr_cols, &nr_nzeros);
    
    ///> Init pthreads.
    tids = (pthread_t *) malloc((A->nr_threads - 1) * sizeof(pthread_t));
    if (!tids) {
	fprintf(stderr, "malloc failed\n");
	exit(1);
    } 
    if (pthread_barrier_init(&barrier, NULL, A->nr_threads)) {
        fprintf(stderr, "pthread_barrier_init");
        exit(1);
    }
    
    ///> Determine number of loops for CG.
    nr_loops = atoi(argv[1]);
    
    ///> Create vectors with the appropriate size.
    x = vector_double_create(nr_cols);
    sol = vector_double_create(nr_cols);
    b = vector_double_create(nr_rows);
    r = vector_double_create(nr_cols);
    p = vector_double_create(nr_cols);
    t = vector_double_create(nr_cols);
    
    ///> Assign the appropriate parameters to each thread.
    cg_params *params = (cg_params *) malloc(A->nr_threads * sizeof(cg_params));
    for (i = 0; i < A->nr_threads; i++) {
        params[i].nr_loops = nr_loops;
        params[i].spm_thread = A->spm_threads+i;
        params[i].in = p;
        params[i].out = t;
    }
    
    ///> Load vector sol with random values.
    vector_double_init_rand_range(sol, (double) -1000, (double) 1000);
	
    ///> Do A*sol=b to define b vector.
    spm_mt_thread_t *spm_thread = A->spm_threads;
    spm_crs32_double_mt_t *csr_mt = (spm_crs32_double_mt_t *) spm_thread->spm;
    spm_crs32_double_multiply((void *) csr_mt->crs, sol, b);
    sol_dis = cblas_ddot(nr_cols, sol->elements, 1, sol->elements, 1);
    sol_dis = sqrt(sol_dis);
    
    vector_double_init(x, 0.01);                                ///> Initiate vector x (x0).
    spm_crs32_double_multiply((void *) csr_mt->crs, x, t);      ///> Do t = A * x0.
    vector_double_sub(b, t, r);                                 ///> Do r0 = b - t.
    vector_double_copy(r, p);                                   ///> Do p0 = r0.

    ///> Start counting time.
    gettimeofday(&start, NULL);
    
    ///> Set thread to the appropriate cpu.
    setaffinity_oncpu(params[0].spm_thread->cpu);
    
    ///> Initiate side threads.
    for (i = 1; i < A->nr_threads; i++)
	pthread_create(&tids[i-1], NULL, cg_side_thread, (void *) &params[i]);

    for (i = 0; i < nr_loops; i++) {
        pthread_barrier_wait(&barrier);
        spm_crs32_double_mt_multiply((void *) &csr_mt->crs, p, t);      ///> Do t = Ap.
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
        vector_double_mul(p, bi, t);                                    ///> Do p = r + bi*p.
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
    //spm_crs32_double_mt_destroy((void *) A);
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4

