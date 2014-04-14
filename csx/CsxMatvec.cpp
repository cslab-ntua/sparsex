/*
 * CsxMatvec.cpp -- Multithreaded kernel y <-- alpha*A*x + beta*y
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "CsxMatvec.hpp"
#include "ThreadPool.hpp"

using namespace boost;

static spx_vector_t **temp;
static double ALPHA, BETA;
static size_t nr_threads;
#ifdef DISABLE_POOL
static barrier *bar;
#endif

void do_matvec_thread(void *args)
{
	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_double_fn_t *fn = (spmv_double_fn_t *) spm_thread->spmv_fn;
    int start = spm_thread->row_start;
    int end = start + spm_thread->nr_rows;

	if (BETA != 1)
        vec_scale_part(spm_thread->y, spm_thread->y, BETA, start, end);
	fn(spm_thread->spm, spm_thread->x, spm_thread->y, ALPHA);
}

void do_matvec_sym_thread(void *args)
{
	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_double_sym_fn_t *fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
    int start = spm_thread->row_start;
    int end = start + spm_thread->nr_rows;
	int id = spm_thread->id;
    bool *sense = spm_thread->sense;

    vec_init_from_map(temp, 0, spm_thread->map);
    central_barrier(sense, nr_threads);
	if (BETA != 1)
        vec_scale_part(spm_thread->y, spm_thread->y, BETA, start, end);
    fn(spm_thread->spm, spm_thread->x, spm_thread->y, temp[id], ALPHA);
    central_barrier(sense, nr_threads);
    /* Switch Reduction Phase */
    vec_add_from_map(spm_thread->y, temp, spm_thread->y, spm_thread->map);
}

void MatVecMult(spm_mt_t *spm_mt, spx_vector_t *x, spx_scalar_t alpha,
                spx_vector_t *y, spx_scalar_t beta)
{
    ALPHA = alpha; BETA = beta;
    nr_threads = spm_mt->nr_threads;

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
	}

#ifdef DISABLE_POOL
    vector<std::shared_ptr<thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<thread>
            (bind(&do_matvec_thread, spm_mt->spm_threads + i));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#else
    ThreadPool &pool = ThreadPool::GetInstance();
    for (size_t i = 0; i < pool.GetSize(); i++) {
        pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
    }
    pool.SetKernel(SPMV);
    central_barrier(pool.GetSense(), nr_threads);
    do_matvec_thread(spm_mt->spm_threads);
    central_barrier(pool.GetSense(), nr_threads);
#endif
}

void MatVecMult_sym(spm_mt_t *spm_mt, spx_vector_t *x, spx_scalar_t alpha, 
                    spx_vector_t *y, spx_scalar_t beta)
{
    ALPHA = alpha; BETA = beta;
	nr_threads = spm_mt->nr_threads;
    temp = spm_mt->local_buffers;
	temp[0] = y;
#ifdef DISABLE_POOL
    bar = new barrier(nr_threads);
#endif

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
	}

#ifdef DISABLE_POOL
    vector<std::shared_ptr<thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<thread>
            (boost::bind(&do_matvec_sym_thread, spm_mt->spm_threads + i));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#else
    ThreadPool &pool = ThreadPool::GetInstance();
    for (size_t i = 0; i < pool.GetSize(); i++) {
        pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
    }
    spm_mt->spm_threads[0].sense = pool.GetSense();
    pool.SetKernel(SPMV_SYM);
    central_barrier(pool.GetSense(), nr_threads);
    do_matvec_sym_thread(spm_mt->spm_threads);
    central_barrier(pool.GetSense(), nr_threads);
#endif
}
