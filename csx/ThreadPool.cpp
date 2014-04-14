/*
 * ThreadPool.cpp --  Thread management in SparseX.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "ThreadPool.hpp"

ThreadPool::~ThreadPool() 
{
    work_done_.store(true);
    central_barrier(GetSense(), size_ + 1);
    central_barrier(GetSense(), size_ + 1);

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

    central_barrier(GetSense(), size_ + 1);
}

void ThreadPool::Run(Worker *worker)
{
    central_barrier(worker->GetSense(), size_ + 1);

    // Wait for the main thread to set a kernel to be executed
    central_barrier(worker->GetSense(), size_ + 1);
    while (!work_done_.load()) {
        switch (worker->GetJob()) {
        case SPMV:
            do_matvec_thread(worker->data_);
            break;
        case SPMV_SYM:
            do_matvec_sym_thread(worker->data_);
            break;
        default:
            break;
        }
        central_barrier(worker->GetSense(), size_ + 1);
        // Wait for a new kernel to be set
        central_barrier(worker->GetSense(), size_ + 1);
    }
    central_barrier(worker->GetSense(), size_ + 1);
}

void ThreadPool::SetKernel(int kernel_id)
{
    for (size_t i = 0; i < size_; i++) {
        workers_[i].SetJob(kernel_id);
        if (kernel_id == SPMV_SYM) {
            spm_mt_thread_t *data = (spm_mt_thread_t *) workers_[i].data_;
            data->sense = workers_[i].GetSense();
        }
    }
}
