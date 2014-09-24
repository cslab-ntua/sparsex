/*
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file CsxKernels.cpp
 * \brief Multithreaded SpMV kernels
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/CsxKernels.hpp>
#include <iostream>
using namespace sparsex;
using namespace sparsex::runtime;

vector_t **temp;
size_t nr_threads;

#if SPX_DISABLE_POOL
using namespace boost;
static barrier *bar;
#endif

void MatVecKernel(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha,
                  vector_t *y, spx_value_t beta)
{
    nr_threads = spm_mt->nr_threads;

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
		spm_mt->spm_threads[i].alpha = alpha;
		spm_mt->spm_threads[i].beta = beta;
	}

#if SPX_DISABLE_POOL
    vector<std::shared_ptr<thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<thread>
            (bind(&do_kernel_thread, spm_mt->spm_threads + i));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#else
    ThreadPool &pool = ThreadPool::GetInstance();
    for (size_t i = 0; i < pool.GetSize(); i++) {
        pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
    }
    pool.SetKernel(SPMV_KERNEL);
    centralized_barrier(pool.GetSense(), nr_threads);
    do_kernel_thread(spm_mt->spm_threads);
    centralized_barrier(pool.GetSense(), nr_threads);
#endif
}

void MatVecKernel_sym(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha,
                      vector_t *y, spx_value_t beta)
{
	nr_threads = spm_mt->nr_threads;
    temp = spm_mt->local_buffers;
	temp[0] = y;
#if SPX_DISABLE_POOL
    bar = new barrier(nr_threads);
#endif

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
		spm_mt->spm_threads[i].alpha = alpha;
		spm_mt->spm_threads[i].beta = beta;
	}

#if SPX_DISABLE_POOL
    vector<std::shared_ptr<thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<thread>
            (boost::bind(&do_kernel_sym_thread, spm_mt->spm_threads + i));
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
    pool.SetKernel(SPMV_KERNEL_SYM);
    centralized_barrier(pool.GetSense(), nr_threads);
    do_kernel_sym_thread(spm_mt->spm_threads);
    centralized_barrier(pool.GetSense(), nr_threads);
#endif
}

void MatVecMult(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha,
                vector_t *y)
{
    nr_threads = spm_mt->nr_threads;

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
		spm_mt->spm_threads[i].alpha = alpha;
	}

    VecInit(y, 0);

#if SPX_DISABLE_POOL
    vector<std::shared_ptr<thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<thread>
            (bind(&do_mv_thread, spm_mt->spm_threads + i));
    }

    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i]->join();
    }
#else
    ThreadPool &pool = ThreadPool::GetInstance();
    for (size_t i = 0; i < pool.GetSize(); i++) {
        pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
    }
    pool.SetKernel(SPMV_MULT);
    centralized_barrier(pool.GetSense(), nr_threads);
    do_mv_thread(spm_mt->spm_threads);
    centralized_barrier(pool.GetSense(), nr_threads);
#endif
}

void MatVecMult_sym(spm_mt_t *spm_mt, vector_t *x, spx_value_t alpha, 
                    vector_t *y)
{
	nr_threads = spm_mt->nr_threads;
    temp = spm_mt->local_buffers;
	temp[0] = y;
#if SPX_DISABLE_POOL
    bar = new barrier(nr_threads);
#endif

	for (size_t i = 0; i < nr_threads; i++) {
		spm_mt->spm_threads[i].x = x;
		spm_mt->spm_threads[i].y = y;
		spm_mt->spm_threads[i].alpha = alpha;
	}

    VecInit(y, 0);

#if SPX_DISABLE_POOL
    vector<std::shared_ptr<thread> > threads(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        threads[i] = make_shared<thread>
            (boost::bind(&do_mv_sym_thread, spm_mt->spm_threads + i));
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
    pool.SetKernel(SPMV_MULT_SYM);
    centralized_barrier(pool.GetSense(), nr_threads);
    do_mv_sym_thread(spm_mt->spm_threads);
    centralized_barrier(pool.GetSense(), nr_threads);
#endif
}
