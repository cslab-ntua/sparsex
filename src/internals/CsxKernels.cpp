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

namespace sparsex {
  namespace runtime {
    extern atomic<int> global_sense;
  }
}

void MatVecKernel(const spm_mt_t *spm_mt, const vector_t *x,
                  spx_value_t alpha, vector_t *y, spx_value_t beta)
{
  nr_threads = spm_mt->nr_threads;

  for (size_t i = 0; i < nr_threads; i++) {
    spm_mt->spm_threads[i].x = x;
    spm_mt->spm_threads[i].y = y;
    spm_mt->spm_threads[i].alpha = alpha;
    spm_mt->spm_threads[i].beta = beta;
  }

  ThreadPool &pool = ThreadPool::GetInstance();
  for (size_t i = 0; i < pool.GetSize(); i++) {
    pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
  }
  pool.SetKernel(SPMV_KERNEL);
  centralized_barrier(pool.GetSense(), nr_threads);
  do_kernel_thread(spm_mt->spm_threads);
  centralized_barrier(pool.GetSense(), nr_threads);
}

void MatVecKernel_sym(const spm_mt_t *spm_mt, const vector_t *x,
                      spx_value_t alpha, vector_t *y, spx_value_t beta)
{
  nr_threads = spm_mt->nr_threads;
  temp = spm_mt->local_buffers;
  temp[0] = y;

  for (size_t i = 0; i < nr_threads; i++) {
    spm_mt->spm_threads[i].x = x;
    spm_mt->spm_threads[i].y = y;
    spm_mt->spm_threads[i].alpha = alpha;
    spm_mt->spm_threads[i].beta = beta;
  }

  ThreadPool &pool = ThreadPool::GetInstance();
  for (size_t i = 0; i < pool.GetSize(); i++) {
    pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
  }
  spm_mt->spm_threads[0].sense = !global_sense;
  pool.SetKernel(SPMV_KERNEL_SYM);
  centralized_barrier(pool.GetSense(), nr_threads);
  do_kernel_sym_thread(spm_mt->spm_threads);
  centralized_barrier(pool.GetSense(), nr_threads);
}

void MatVecMult(const spm_mt_t *spm_mt, const vector_t *x,
		spx_value_t alpha, vector_t *y)
{
  nr_threads = spm_mt->nr_threads;

  for (size_t i = 0; i < nr_threads; i++) {
    spm_mt->spm_threads[i].x = x;
    spm_mt->spm_threads[i].y = y;
    spm_mt->spm_threads[i].alpha = alpha;
  }

  VecInit(y, 0);

  ThreadPool &pool = ThreadPool::GetInstance();
  for (size_t i = 0; i < pool.GetSize(); i++) {
    pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
  }
  pool.SetKernel(SPMV_MULT);
  centralized_barrier(pool.GetSense(), nr_threads);
  do_mv_thread(spm_mt->spm_threads);
  centralized_barrier(pool.GetSense(), nr_threads);
}

void MatVecMult_sym(const spm_mt_t *spm_mt, const vector_t *x,
                    spx_value_t alpha, vector_t *y)
{
  nr_threads = spm_mt->nr_threads;
  temp = spm_mt->local_buffers;
  temp[0] = y;

  for (size_t i = 0; i < nr_threads; i++) {
    spm_mt->spm_threads[i].x = x;
    spm_mt->spm_threads[i].y = y;
    spm_mt->spm_threads[i].alpha = alpha;
  }

  VecInit(y, 0);

  ThreadPool &pool = ThreadPool::GetInstance();
  for (size_t i = 0; i < pool.GetSize(); i++) {
    pool.SetWorkerData(i, spm_mt->spm_threads + i + 1);
  }
  spm_mt->spm_threads[0].sense = !global_sense;
  pool.SetKernel(SPMV_MULT_SYM);
  centralized_barrier(pool.GetSense(), nr_threads);
  do_mv_sym_thread(spm_mt->spm_threads);
  centralized_barrier(pool.GetSense(), nr_threads);
}
