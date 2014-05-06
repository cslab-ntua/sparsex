/*
 * ThreadPool.cpp --  Thread management in SparseX.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "sparsex/internals/ThreadPool.hpp"

ThreadPool::~ThreadPool() 
{
    work_done_.store(true);
    centralized_barrier(GetSense(), size_ + 1);
    centralized_barrier(GetSense(), size_ + 1);

    for (size_t i = 0; i < size_; i++) {
        if (workers_[i].thread_->joinable())
            workers_[i].thread_->join();
    }
}

void ThreadPool::InitThreads(size_t nr_threads)
{
    if (nr_threads < 0) {
        LOG_ERROR << "initializing thread pool with unvalid number of threads";
        exit(1);
    }

    size_ = nr_threads;
    barrier_cnt = nr_threads + 1;
    if (size_)
        workers_.resize(size_);

    for (size_t i = 0; i < size_; i++) {
        workers_[i].SetId(i);
        // You can pass a member function pointer as the function, provided you
        // supply a suitable object pointer as the first argument
        workers_[i].thread_ = make_shared<boost::thread>
            (&ThreadPool::Run, this, workers_.data() + i);  // exception throwing try-catch
    }

    centralized_barrier(GetSense(), size_ + 1);
}

void ThreadPool::Run(Worker *worker)
{
    centralized_barrier(worker->GetSense(), size_ + 1);

    // Wait for the main thread to set a kernel to be executed
    centralized_barrier(worker->GetSense(), size_ + 1);
    while (!work_done_.load()) {
        switch (worker->GetJob()) {
        case SPMV_MULT:
            do_mv_thread(worker->data_);
            break;
        case SPMV_KERNEL:
            do_kernel_thread(worker->data_);
            break;
        case SPMV_MULT_SYM:
            do_mv_sym_thread(worker->data_);
            break;
        case SPMV_KERNEL_SYM:
            do_kernel_sym_thread(worker->data_);
            break;
        default:
            break;
        }
        centralized_barrier(worker->GetSense(), size_ + 1);
        // Wait for a new kernel to be set
        centralized_barrier(worker->GetSense(), size_ + 1);
    }
    centralized_barrier(worker->GetSense(), size_ + 1);
}

void ThreadPool::SetKernel(int kernel_id)
{
    for (size_t i = 0; i < size_; i++) {
        workers_[i].SetJob(kernel_id);
        if (kernel_id == SPMV_MULT_SYM || kernel_id == SPMV_KERNEL_SYM) {
            spm_mt_thread_t *data = (spm_mt_thread_t *) workers_[i].data_;
            data->sense = workers_[i].GetSense();
        } else if (kernel_id == IDLE) {
            workers_[i].data_ = NULL;
        }
    }
}
