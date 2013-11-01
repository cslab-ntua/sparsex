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
	if (BETA != 1)
        vec_scale_part(spm_thread->y, spm_thread->y, BETA, start, end);
	fn(spm_thread->spm, spm_thread->x, spm_thread->y, ALPHA);
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
	if (BETA != 1)
        vec_scale_part(spm_thread->y, spm_thread->y, BETA, start, end);
    // fn(spm_thread->spm, spm_thread->x, params->tmp, temp[id]);
    fn(spm_thread->spm, spm_thread->x, spm_thread->y, temp[id]);
    cur_barrier.wait();
    /* Switch Reduction Phase */
    // vec_add_from_map(params->tmp, temp, params->tmp, spm_thread->map);
    vec_add_from_map(spm_thread->y, temp, spm_thread->x, spm_thread->map);
    cur_barrier.wait();
    // if (ALPHA != 1)
    //     vec_scale_part(params->tmp, params->tmp, ALPHA, start, end);
    // vec_scale_add_part(params->tmp, spm_thread->y, spm_thread->y, BETA,
    //                    start, end);
    // cur_barrier.wait();
}

void matvec_mt(spm_mt_t *spm_mt, vector_t *x, double alpha, vector_t *y,
               double beta)
{
    ALPHA = alpha; BETA = beta;
	nr_threads = spm_mt->nr_threads;
    boost::barrier cur_barrier(nr_threads);
    matvec_params *params = (matvec_params *) malloc(sizeof(matvec_params) *
                                                     nr_threads);
    size_t rows = 0;
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
        params[i].spm_thread = spm_mt->spm_threads + i;
        params[i].start = rows;
        params[i].end = rows + spm_mt->spm_threads[i].nr_rows;
        rows += spm_mt->spm_threads[i].nr_rows;
	}

    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
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
		int node = spm_mt->spm_threads[i].node;
		temp[i] = vec_create_onnode(y->size, node);
	}
#else
	for (unsigned int i = 1; i < nr_threads; i++)
		temp[i] = vec_create(n, NULL);
#endif

    size_t rows = 0;
	for (unsigned int i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
        params[i].spm_thread = spm_mt->spm_threads + i;
        params[i].start = rows;
        params[i].end = rows + spm_mt->spm_threads[i].nr_rows;
        rows += spm_mt->spm_threads[i].nr_rows;
	}

    std::vector<boost::shared_ptr<boost::thread> > threads(nr_threads);
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
    free(params);
}
