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

static spx_vector_t **temp;
static double ALPHA = 1, BETA = 1;

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

void do_matvec_sym_thread(void *args, boost::barrier &cur_barrier)
{
	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) args;
	setaffinity_oncpu(spm_thread->cpu);
	spmv_double_sym_fn_t *fn = (spmv_double_sym_fn_t *) spm_thread->spmv_fn;
    int start = spm_thread->row_start;
    int end = start + spm_thread->nr_rows;
	int id = spm_thread->id;

    vec_init_from_map(temp, 0, spm_thread->map);
    cur_barrier.wait();
	if (BETA != 1)
        vec_scale_part(spm_thread->y, spm_thread->y, BETA, start, end);
    fn(spm_thread->spm, spm_thread->x, spm_thread->y, temp[id], ALPHA);
    cur_barrier.wait();
    /* Switch Reduction Phase */
    vec_add_from_map(spm_thread->y, temp, spm_thread->y, spm_thread->map);
}

void MatVecMult(spm_mt_t *spm_mt, spx_vector_t *x, spx_scalar_t alpha,
                spx_vector_t *y, spx_scalar_t beta)
{
    ALPHA = alpha; BETA = beta;
	size_t nr_threads = spm_mt->nr_threads;

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
	}

#ifdef USE_THREAD_POOL
    ThreadPool &pool = ThreadPool::GetInstance();
    for (size_t i = 0; i < pool.GetSize(); i++) {
        pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
    }
    pool.SetKernel(SPMV);
#ifdef USE_BOOST_BARRIER
    pool.start_barrier_->wait();
#else
    central_barrier(pool.local_sense, nr_threads);
#endif

    // tsc_start(&thread_tsc);
    do_matvec_thread(spm_mt->spm_threads);
    // tsc_pause(&thread_tsc);

#ifdef USE_BOOST_BARRIER
    pool.end_barrier_->wait();
#else
    central_barrier(pool.local_sense, nr_threads);
#endif

#else
    vector<shared_ptr<boost::thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<boost::thread>
            (boost::bind(&do_matvec_thread, spm_mt->spm_threads + i));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#endif
}

void MatVecMult_sym(spm_mt_t *spm_mt, spx_vector_t *x, spx_scalar_t alpha, 
                    spx_vector_t *y, spx_scalar_t beta)
{
    ALPHA = alpha; BETA = beta;
	size_t nr_threads = spm_mt->nr_threads;
    temp = spm_mt->local_buffers;
	temp[0] = y;
    boost::barrier cur_barrier(nr_threads);

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
	}

    vector<shared_ptr<boost::thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<boost::thread>
            (boost::bind(&do_matvec_sym_thread, spm_mt->spm_threads + i,
                         boost::ref(cur_barrier)));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
}
