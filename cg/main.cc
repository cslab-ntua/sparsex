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
#include <unistd.h>

#include <cmath>
#include <cstdlib>

#include "spmv.h"
#include "cg.h"

using namespace std;

///> Type of method used in CG.
typedef enum {
   CSR_SPMV = 0,
   SSS_SPMV,
   CSX_SPMV,
   CSX_SYM_SPMV
} cg_method_t;


int main(int argc, char *argv[])
{
    unsigned long       i;
    long                j;
    char                *mmf_file;
    char                *symmetric;
    
    uint64_t            nrows, ncols, nnz, n, ncpus;
    spm_mt_t            *spm_mt;
    vector_double_t     *x, *b, *sol, *r, *p, *t;
    vector_double_t     **sub_p = NULL;
    double              rd, sol_dis;
    
    double              *rr, *tp, *rr_new;
    double              *ai, *bi;
    
    pthread_t           *tids;
    pthread_barrier_t   barrier;
    
    method_t            *meth;
    
    csx_double_t	*csx;
    csx_double_sym_t    *csx_sym;

    double              cg_time;
    double              spmv_time = 0;
    double              red_time = 0;
    
    ///> Function variables.   
    void                *(*CgSideThread)(void *);
    void                (*CgMainThread)(cg_params *, double *, double *, double *);
    uint64_t            (*SpMSize)(void *);
    
    ///> Default options.
    unsigned long       nLoops = 1;
    unsigned long       nloops = 512;
    cg_method_t         cg_method = CSR_SPMV;
    
    csx::CsxExecutionEngine &engine = csx::CsxJitInit();

    ///> Parse Options.
    while ((j = getopt(argc, argv, "xsl:L:")) != -1) {
        switch (j) {
            case 'x':
                if (cg_method == CSR_SPMV)
                    cg_method = CSX_SPMV;
                else
                    cg_method = CSX_SYM_SPMV;
                break;
            case 's':
                if (cg_method == CSR_SPMV)
                    cg_method = SSS_SPMV;
                else
                    cg_method = CSX_SYM_SPMV;
                break;
            case 'l':
                nloops = atol(optarg);
                break;
            case 'L':
                nLoops = atol(optarg);
                break;
	    default:
	        fprintf(stderr, "Usage: cg [-x -s -l <number of inside loops> -L <number of outside loops>] mmf_file\n");
	        exit(1);
	}
    }
    
    ///> Take input file from which matrix will be loaded.
    if (argc != optind+1) {
        //fprintf(stderr, "Usage: cg [-x -s -l <number of loops> -L <number of outside loops>] mmf_file\n");
        exit(1);
    }
    mmf_file = argv[optind];
    
    ///> Load matrix in appropriate format.
    switch(cg_method) {
        case(CSR_SPMV):
#ifndef SPM_NUMA
            spm_mt = (spm_mt_t *) spm_crs32_double_mt_init_mmf(mmf_file, &nrows,
                                                               &ncols, &nnz,
                                                               NULL);
	    meth = method_get((char *) "spm_crs32_double_mt_multiply");
#else
            spm_mt = (spm_mt_t *) spm_crs32_double_mt_numa_init_mmf(mmf_file,
                                                                    &nrows,
                                                                    &ncols,
                                                                    &nnz, NULL);
	    meth = method_get((char *) "spm_crs32_double_mt_numa_multiply");
#endif
            for (i = 0; i < spm_mt->nr_threads; i++)
                spm_mt->spm_threads[i].spmv_fn = meth->fn;
            
            assert(nrows == ncols && "Matrix is not square");
	    n = nrows;
            
            CgSideThread = NormalCgSideThread;
            CgMainThread = NormalCgMainThread;
            SpMSize = spm_crs32_double_mt_size;
            
            break;
        
        case(SSS_SPMV):
#ifndef SPM_NUMA
            spm_mt = (spm_mt_t *)
                spm_crs32_double_sym_mt_init_mmf(mmf_file, &nrows, &ncols,
                                                 &nnz, NULL);
	    meth = method_get((char *) "spm_crs32_double_sym_mt_multiply");
#else
            spm_mt = (spm_mt_t *)
                spm_crs32_double_sym_mt_numa_init_mmf(mmf_file, &nrows, &ncols,
                                                      &nnz, NULL);
	    meth = method_get((char *) "spm_crs32_double_sym_mt_numa_multiply");
#endif
            for (i = 0; i < spm_mt->nr_threads; i++)
                spm_mt->spm_threads[i].spmv_fn = meth->fn;
            
            assert(nrows == ncols && "Matrix is not square");
            n = nrows;
            
            CgSideThread = SymCgSideThread;
            CgMainThread = SymCgMainThread;
            SpMSize = spm_crs32_double_sym_mt_size;
            
            break;
        
        case(CSX_SPMV) :
            symmetric = getenv("SYMMETRIC");
            assert(symmetric == NULL && 
                   "environment variable SYMMETRIC must not be set");

            spm_mt = GetSpmMt(mmf_file, engine);
            
            csx = (csx_double_t *) spm_mt->spm_threads[0].spm;
            ncols = csx->ncols;
	    nrows = 0;
            nnz = 0;
            for (i = 0; i < spm_mt->nr_threads; i++) {
            	csx = (csx_double_t *) spm_mt->spm_threads[i].spm;
            	nrows += csx->nrows;
		nnz += csx->nnz;
	    }
	    assert(nrows == ncols && "Matrix is not square");
	    n = nrows;
	    
            CgSideThread = NormalCgSideThread;
            CgMainThread = NormalCgMainThread;
            SpMSize = CsxSize;
            
            break;
        
        case(CSX_SYM_SPMV) :
            symmetric = getenv("SYMMETRIC");
            assert(symmetric != NULL && 
                   "environment variable SYMMETRIC must be set");
            
            spm_mt = GetSpmMt(mmf_file, engine);
            
            csx_sym = (csx_double_sym_t *) spm_mt->spm_threads[0].spm;
            csx = (csx_double_t *) csx_sym->lower_matrix;
            ncols = csx->ncols;
	    nrows = 0;
            nnz = 0;
            for (i = 0; i < spm_mt->nr_threads; i++) {
            	csx_sym = (csx_double_sym_t *) spm_mt->spm_threads[i].spm;
                csx = (csx_double_t *) csx_sym->lower_matrix;
            	nrows += csx->nrows;
		nnz += csx->nnz;
	    }
	    assert(nrows == ncols && "Matrix is not square");
	    n = nrows;
	    
            CgSideThread = SymCgSideThread;
            CgMainThread = SymCgMainThread;
            SpMSize = CsxSymSize;
            
            break;
            
        default :
            fprintf(stderr, "Wrong method choosed\n");
            exit(1);
    }
    
    ncpus = spm_mt->nr_threads;

#ifdef SPM_NUMA
    size_t *parts = (size_t *) malloc(ncpus * sizeof(*parts));
    size_t *matrix_parts = (size_t *) malloc(ncpus * sizeof(*matrix_parts));
    int *nodes = (int *) malloc(ncpus * sizeof(*nodes));
    
    for (i = 0; i < ncpus; i++) {
	spm_mt_thread_t *spm_thread = spm_mt->spm_threads + i;
	
	parts[i] = ((i + 1) * n / ncpus - i * n / ncpus) * sizeof(double);
	matrix_parts[i] = spm_thread->nr_rows * sizeof(double);
	nodes[i] = numa_node_of_cpu(spm_thread->cpu);
    }
#endif

    ///> Init pthreads.
    tids = (pthread_t *) malloc((ncpus - 1) * sizeof(pthread_t));
    if (!tids) {
	fprintf(stderr, "Malloc of pthreads failed\n");
	exit(1);
    }
    if (pthread_barrier_init(&barrier, NULL, ncpus)) {
        fprintf(stderr, "Pthread barrier init failed");
        exit(1);
    }

    for (unsigned int loops=0; loops < nLoops; loops++) {  
        ///> Create vectors with the appropriate size.
        sol = vector_double_create(n);
        b = vector_double_create(n);
#ifndef SPM_NUMA
        x = vector_double_create(n);
        r = vector_double_create(n);
        p = vector_double_create(n);
        t = vector_double_create(n);

        if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV) {
            sub_p = (vector_double_t **)
                        malloc(ncpus * sizeof(vector_double_t *));
            sub_p[0] = t;
            for (i = 1; i < ncpus; i++)
                sub_p[i] = vector_double_create(n);
        }
#else
        x = vector_double_create_interleaved(n, parts, ncpus, nodes);
        r = vector_double_create_interleaved(n, parts, ncpus, nodes);
        p = vector_double_create_interleaved(n, parts, ncpus, nodes);
        t = vector_double_create_interleaved(n, matrix_parts, ncpus, nodes);
        
        if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV) {
            sub_p = (vector_double_t **) 
                       malloc(ncpus * sizeof(vector_double_t *));
            sub_p[0] = t;
            for (i = 1; i < ncpus; i++)
                sub_p[i] = vector_double_create_onnode(n, nodes[i]);
        }
#endif
        if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV)
            for (i = 1; i < ncpus; i++) 
                vector_double_init(sub_p[i], 0);
    
        ///> Initialize partial values of doubles.
        rr = (double *) malloc(ncpus * sizeof(double));
        tp = (double *) malloc(ncpus * sizeof(double));
        rr_new = (double *) malloc(ncpus * sizeof(double));
        ai = (double *) malloc(sizeof(double));
        bi = (double *) malloc(sizeof(double));

        ///> Assign the appropriate parameters to each thread.
        cg_params *params = (cg_params *) malloc(ncpus * sizeof(cg_params));
                            
        for (i = 0; i < ncpus; i++) {
            params[i].nloops = nloops;
            params[i].ncpus = ncpus;
            params[i].spm_thread = spm_mt->spm_threads+i;
            params[i].p = p;
            params[i].temp = t;
            params[i].r = r;
            params[i].x = x;
            params[i].barrier = &barrier;
            params[i].start = i * n / ncpus;
            params[i].end = (i + 1) * n / ncpus;
            params[i].rr = rr + i;
            params[i].tp = tp + i;
            params[i].rr_new = rr_new + i;
            params[i].ai = ai;
            params[i].bi = bi;
            if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV)
                params[i].sub_vectors = sub_p;
        }

        ///> Load vector solution with random values and calculate its distance.
        vector_double_init_rand_range(sol, (double) -1000, (double) 1000);
        sol_dis = vector_double_mul(sol, sol);
        sol_dis = sqrt(sol_dis);

        ///> Find b vector for the specified solution.
        if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV)
            FindSymSolution(spm_mt, sol, b, t);
        else
            FindSolution(spm_mt, sol, b, t);
    
        ///> Initialize CG method.
        if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV)
            InitializeSymCg(spm_mt, x, r, p, b, t);
        else
            InitializeCg(spm_mt, x, r, p, b, t);

        ///> Initiate side threads.
        for (i = 1; i < ncpus; i++)
	    pthread_create(&tids[i-1], NULL, CgSideThread, (void *) &params[i]);

        ///> Execute main thread.
        CgMainThread(params, &cg_time, &spmv_time, &red_time);

#ifdef SPM_NUMA	
        int alloc_err;

	alloc_err = check_interleaved(x->elements, parts, ncpus, nodes);
	print_alloc_status("x", alloc_err);

	alloc_err = check_interleaved(r->elements, parts, ncpus, nodes);
	print_alloc_status("r", alloc_err);

	alloc_err = check_interleaved(p->elements, parts, ncpus, nodes);
	print_alloc_status("p", alloc_err);

	alloc_err = check_interleaved(t->elements, matrix_parts, ncpus, nodes);
	print_alloc_status("t", alloc_err);
	
	if (cg_method == SSS_SPMV || cg_method == CSX_SYM_SPMV) {
	    alloc_err = 0;
	    for (i = 1; i < ncpus; i++)
	        alloc_err += check_region(sub_p[i]->elements,
		                          n * sizeof(double), nodes[i]);
	    print_alloc_status("temporary buffers", alloc_err);
        }
#endif

        ///> Find relative distance.
        vector_double_sub(sol, x, t);
        rd = vector_double_mul(t, t);
        rd = sqrt(rd);
        rd /= sol_dis;
    
        ///> Print Results.
        printf("m:%s l:%lu rd:%lf st:%lf rt:%lf ct:%lf\n", basename(mmf_file), 
               nloops, rd, spmv_time, red_time, cg_time);
    
        ///> Release vectors.
        vector_double_destroy(x);
        vector_double_destroy(sol);
        vector_double_destroy(b);
        vector_double_destroy(r);
        vector_double_destroy(p);
        vector_double_destroy(t);
        /*
        if (cg_method == CSR_sym || cg_method == CSX_sym) {
            for (i = 1; i < ncpus; i++)
                vector_double_destroy(sub_p[i]);
            free(sub_p);
        }
        */
        
        ///> Free parameters.
        free(rr);
        free(tp);
        free(rr_new);
        free(ai);
        free(bi);
        free(params);
    
    }
    
    ///> Free pthreads.
    free(tids);
    
    ///> Release matrix.
    free(spm_mt->spm_threads);
    free(spm_mt);
    
    return 0;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
