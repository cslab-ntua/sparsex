/* -*- C++ -*-
 *
 * CsxMatvec.cc -- Multithreaded kernel y <-- alpha*A*x + beta*y
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "CsxMatvec.hpp"

static unsigned int nr_threads = 0;
static vector_t **temp = NULL;
static double ALPHA = 1, BETA = 1;

static void do_matvec_thread(matvec_params *params, boost::barrier &cur_barrier)
{
	spm_mt_thread_t *spm_thread = params->spm_thread;
	spmv_double_fn_t *fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
    int start = params->start;
    int end = params->end;

	setaffinity_oncpu(spm_thread->cpu);

    cur_barrier.wait();
	fn(spm_thread->spm, spm_thread->x, params->tmp);
	if (ALPHA != 1)
        vec_scale_part(params->tmp, params->tmp, ALPHA, start, end);
	vec_scale_add_part(params->tmp, spm_thread->y, spm_thread->y, BETA,
                       start, end);
    cur_barrier.wait();
}

static void do_matvec_sym_thread(matvec_params *params,
                                 boost::barrier &cur_barrier)
{
	spm_mt_thread_t *spm_thread = params->spm_thread;
	spmv_double_sym_fn_t *fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
    int start = params->start;
    int end = params->end;
	int id = spm_thread->id;

	setaffinity_oncpu(spm_thread->cpu);

    vec_init_from_map(temp, 0, spm_thread->map);
    cur_barrier.wait();
    fn(spm_thread->spm, spm_thread->x, params->tmp, temp[id]);
    cur_barrier.wait();
    /* Switch Reduction Phase */
    vec_add_from_map(params->tmp, temp, params->tmp, spm_thread->map);
    cur_barrier.wait();
    if (ALPHA != 1)
        vec_scale_part(params->tmp, params->tmp, ALPHA, start, end);
    vec_scale_add_part(params->tmp, spm_thread->y, spm_thread->y, BETA,
                       start, end);
    cur_barrier.wait();
}

void matvec_mt(spm_mt_t *spm_mt, vector_t *x, double alpha, vector_t *y,
               double beta)
{
    ALPHA = alpha; BETA = beta;
	nr_threads = spm_mt->nr_threads;
    boost::barrier cur_barrier(nr_threads);
    matvec_params *params = (matvec_params *) malloc(sizeof(matvec_params) *
                                                     nr_threads);

#ifdef SPM_NUMA
	size_t *xparts = (size_t *) malloc(sizeof(*xparts) * nr_threads);
	size_t *yparts = (size_t *) malloc(sizeof(*yparts) * nr_threads);
	size_t *tparts = (size_t *) malloc(sizeof(*tparts) * nr_threads);
	int *xnodes = (int *) malloc(sizeof(*xnodes) * nr_threads);
	int *ynodes = (int *) malloc(sizeof(*ynodes) * nr_threads);
	int *tnodes = (int *) malloc(sizeof(*tnodes) * nr_threads);
	if (!xparts || !xnodes) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;
		xparts[i] = spm->nr_rows * sizeof(double);
		yparts[i] = spm->nr_rows * sizeof(double);
		tparts[i] = spm->nr_rows * sizeof(double);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
		tnodes[i] = spm->node;
	}

	/* Allocate an interleaved x */
	int alloc_err = 0;
	vector_t *x_int = vec_create_interleaved(x->size, xparts, nr_threads, xnodes);
	vec_copy(x, x_int);
	alloc_err = check_interleaved(x_int->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);

	/* Allocate an interleaved y */
	vector_t *y_int = vec_create_interleaved(y->size, yparts, nr_threads, ynodes);
	vec_copy(y, y_int);
	alloc_err = check_interleaved(y_int->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);

	/* Allocate an interleaved tmp */
	vector_t *tmp = vec_create_interleaved(y->size, tparts, nr_threads, tnodes);
	vec_init(tmp, 0);
	alloc_err = check_interleaved(tmp->elements, tparts, nr_threads, tnodes);
	print_alloc_status("temp vector", alloc_err);

	size_t rows = 0;
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x_int;
		spm_mt->spm_threads[i].y = y_int;
        params[i].spm_thread = spm_mt->spm_threads + i;
        params[i].tmp = tmp;
        params[i].start = rows;
        params[i].end = rows + spm_mt->spm_threads[i].nr_rows;
        rows += spm_mt->spm_threads[i].nr_rows;
	}
#else
    size_t rows = 0;
    vector_t *tmp = vec_create(y->size);
	vec_init(tmp, 0);
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
        params[i].spm_thread = spm_mt->spm_threads + i;
        params[i].tmp = tmp;
        params[i].start = rows;
        params[i].end = rows + spm_mt->spm_threads[i].nr_rows;
        rows += spm_mt->spm_threads[i].nr_rows;
	}
#endif

    vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (boost::bind(&do_matvec_thread,
                         params + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

    /* Cleanup */
#ifdef SPM_NUMA
    vec_copy(y_int, y);
    vec_destroy(x_int);
    vec_destroy(y_int);
	free(xparts);
	free(xnodes);
	free(yparts);
	free(ynodes);
	free(tparts);
	free(tnodes);
#endif
    vec_destroy(tmp);
    free(params);
}

void matvec_sym_mt(spm_mt_t *spm_mt, vector_t *x, double alpha, vector_t *y,
                   double beta)
{
    ALPHA = alpha; BETA = beta;
	nr_threads = spm_mt->nr_threads;
    boost::barrier cur_barrier(nr_threads);
	unsigned long n = y->size;
    matvec_params *params = (matvec_params *) malloc(sizeof(matvec_params) *
                                                     nr_threads);

	temp = (vector_t **) malloc(nr_threads * sizeof(vector_t *));
	if (!temp) {
		perror("malloc");
		exit(1);
	}
	temp[0] = y;
#ifdef SPM_NUMA
	for (unsigned int i = 1; i < nr_threads; i++) {
		int tnode = spm_mt->spm_threads[i].node;
		temp[i] = vec_create_onnode(y->size, tnode);
	}
#else
	for (unsigned int i = 1; i < nr_threads; i++)
		temp[i] = vec_create(n);
#endif

#ifdef SPM_NUMA
	size_t *xparts = (size_t *) malloc(sizeof(*xparts) * nr_threads);
	size_t *yparts = (size_t *) malloc(sizeof(*yparts) * nr_threads);
	size_t *tparts = (size_t *) malloc(sizeof(*tparts) * nr_threads);
	int *xnodes = (int *) malloc(sizeof(*xnodes) * nr_threads);
	int *ynodes = (int *) malloc(sizeof(*ynodes) * nr_threads);
	int *tnodes = (int *) malloc(sizeof(*tnodes) * nr_threads);
	if (!xparts || !xnodes) {
		perror("malloc");
		exit(1);
	}

	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt_thread_t *spm = spm_mt->spm_threads + i;
		xparts[i] = spm->nr_rows * sizeof(double);
		yparts[i] = spm->nr_rows * sizeof(double);
		tparts[i] = spm->nr_rows * sizeof(double);
		xnodes[i] = spm->node;
		ynodes[i] = spm->node;
		tnodes[i] = spm->node;
	}

	/* Allocate an interleaved x */
	int alloc_err = 0;
	vector_t *x_int = vec_create_interleaved(n, xparts, nr_threads, xnodes);
	vec_copy(x, x_int);
	alloc_err = check_interleaved(x_int->elements, xparts, nr_threads, xnodes);
	print_alloc_status("input vector", alloc_err);

	/* Allocate an interleaved y */
	vector_t *y_int = vec_create_interleaved(n, yparts, nr_threads, ynodes);
	vec_copy(y, y_int);
	alloc_err = check_interleaved(y_int->elements, yparts, nr_threads, ynodes);
	print_alloc_status("output vector", alloc_err);

	/* Allocate an interleaved tmp */
	vector_t *tmp = vec_create_interleaved(n, tparts, nr_threads, tnodes);
	vec_init(tmp, 0);
	alloc_err = check_interleaved(tmp->elements, tparts, nr_threads, tnodes);
	print_alloc_status("temp vector", alloc_err);

    size_t rows = 0;
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x_int;
		spm_mt->spm_threads[i].y = y_int;
        params[i].spm_thread = spm_mt->spm_threads + i;
        params[i].tmp = tmp;
        params[i].start = rows;
        params[i].end = rows + spm_mt->spm_threads[i].nr_rows;
        rows += spm_mt->spm_threads[i].nr_rows;
	}
#else
    size_t rows = 0;
    vector_t *tmp = vec_create(n);
	vec_init(tmp, 0);
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
        params[i].spm_thread = spm_mt->spm_threads + i;
        params[i].tmp = tmp;
        params[i].start = rows;
        params[i].end = rows + spm_mt->spm_threads[i].nr_rows;
        rows += spm_mt->spm_threads[i].nr_rows;
	}
#endif

    vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = boost::make_shared<boost::thread>
            (boost::bind(&do_matvec_sym_thread,
                         params + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }

    /* Cleanup */
	for (unsigned int i = 1; i < nr_threads; i++)
		vec_destroy(temp[i]);
	free(temp);
    vec_destroy(tmp);
    free(params);
#ifdef SPM_NUMA
    vec_copy(y_int, y);
    vec_destroy(x_int);
    vec_destroy(y_int);
	free(xparts);
	free(xnodes);
	free(yparts);
	free(ynodes);
	free(tparts);
	free(tnodes);
#endif
}
