/*
 * main.cc -- Implementation of Conjugate Gradient Method (CG).
 *
 * Copyright (C) 2011,      Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <pthread.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>

#include <cmath>
#include <cstdlib>

#include "spmv.h"

#include "cg_vector.h"
#include "cg_csr_sym.h"
#include "cg.h"

using namespace std;

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
    int			j;
    char                *mmf_file;
    uint64_t            nrows, ncols, nnz, n;
    spm_mt_t            *spm_mt;
    spm_mt_thread_t     *spm_thread;
    vector_double_t     *x, *b, *sol, *r, *p, *t;
    double              rr, tp, ai, bi, rr_new;
    double              time, rd, sol_dis;
    pthread_t           *tids;
    pthread_barrier_t   barrier;
    spmv_double_fn_t    *fn;
    method_t            *meth;
    csx_double_t	*csx;
    struct timeval      start, end;
    spm_crs32_double_sym_mt_t *crs_mt;
    
    ///> Default options.
    unsigned long   nloops = 512;
    cg_method_t     cg_method = CSR;
    
    ///> Parse Options.
    while ((j = getopt(argc, argv, "cxsl:")) != -1) {
        switch (j) {
            case 'x':
                cg_method = CSX;
                break;
            case 's':
                if (cg_method == CSR)
                    cg_method = CSR_s;
                else
                    cg_method = CSX_s;
                break;
            case 'l':
                nloops = atol(optarg);
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
            spm_mt = (spm_mt_t *) spm_crs32_double_mt_init_mmf(mmf_file, &nrows, &ncols, &nnz);
	    assert(nrows == ncols && "Matrix is not square");
	    n = nrows;
	    meth = method_get((char *) "spm_crs32_double_mt_multiply");
            for (i = 0; i < spm_mt->nr_threads; i++)
                spm_mt->spm_threads[i].spmv_fn = meth->fn;
            break;
        case(CSR_s) :
            spm_mt = (spm_mt_t *) spm_crs32_double_sym_mt_init_mmf(mmf_file, &n, &nnz);
            for (i = 0; i < spm_mt->nr_threads; i++) {
                printf("Thread %lu\n", i);
                crs_mt = (spm_crs32_double_sym_mt_t *) spm_mt->spm_threads[i].spm;
                printf("Start Point (%u, %u)\n", crs_mt->row_start, crs_mt->elem_start);
                printf("End Point (%u, %u)\n", crs_mt->row_end, crs_mt->elem_end);
            }
            exit(0);
            break;
        case(CSX) :
            spm_mt = GetSpmMt(mmf_file);
            csx = (csx_double_t *) spm_mt->spm_threads[0].spm;
            ncols = csx->ncols;
	    nrows = 0;
            nnz = 0;
            for (i = 0; i < spm_mt->nr_threads; i++) {
            	csx = (csx_double_t *) spm_mt->spm_threads[i].spm;
            	nrows += csx->nrows;
		nnz += csx->nnz;
	    }
	    assert(nrows == ncols && "Matrix is not symmetric");
	    n = nrows;
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
    x = vector_double_create(n);
    sol = vector_double_create(n);
    b = vector_double_create(n);
    r = vector_double_create(n);
    p = vector_double_create(n);
    t = vector_double_create(n);

    ///> Assign the appropriate parameters to each thread.
    cg_params *params = (cg_params *) malloc(spm_mt->nr_threads * sizeof(cg_params));
    for (i = 0; i < spm_mt->nr_threads; i++) {
        params[i].nloops = nloops;
        params[i].spm_thread = spm_mt->spm_threads+i;
        params[i].in = p;
        params[i].out = t;
        params[i].barrier = &barrier;
    }
    
    ///> Load vector sol with random values.
    vector_double_init_rand_range(sol, (double) -1000, (double) 1000);
	
    ///> Do A*sol=b to define b vector.
    vector_double_init(b, 0);
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
        vector_double_init(t, 0);
        fn(spm_thread->spm, sol, t);
        vector_double_add(b, t, b);
    }
    sol_dis = vector_double_mul(sol, sol);
    sol_dis = sqrt(sol_dis);
    
    vector_double_init(x, 0.01);                                ///> Initiate vector x (x0).
    vector_double_init(t, 0);                                   ///> Do t = A * x0.
    for (i = 0; i < spm_mt->nr_threads; i++) {
        spm_thread = spm_mt->spm_threads + i;
        fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
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
    fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
    for (i = 0; i < nloops; i++) {
        vector_double_init(t, 0);
        pthread_barrier_wait(&barrier);
        fn(spm_thread->spm, p, t);                                      ///> Do t = Ap.
        pthread_barrier_wait(&barrier);
        rr = vector_double_mul(r, r); 					//Calculate ai.
        tp = vector_double_mul(t, p);
        ai = rr / tp;
        vector_double_scale(t, ai, t);                                  ///> Do r = r - ai*A*p.
        vector_double_sub(r, t, r);
        vector_double_scale(p, ai, t);       				///> Do x = x + ai*p.
	vector_double_add(x, t, x);
	rr_new = vector_double_mul(r, r);   				///> Calculate bi.
	bi = rr_new / rr;
	vector_double_scale(p, bi, t);                                  ///> Do p = r + bi*p.
        vector_double_add(r, t, p);
    }
    
    ///> Stop counting time.
    gettimeofday(&end, NULL);
    time = end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0;
    
    ///> Print results.
    vector_double_sub(sol, x, t);
    rd = vector_double_mul(t, t);
    rd = sqrt(rd);
    rd /= sol_dis;
    printf("m:%s l:%lu rd:%lf t:%lf\n", basename(mmf_file), nloops, rd, time);
    
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
    free(spm_mt->spm_threads);
    free(spm_mt);
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
