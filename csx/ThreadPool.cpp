/*
 * ThreadPool.hpp --  Thread management in SparseX.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "ThreadPool.hpp"

atomic<int> cnt;
atomic<bool> sense(true);
static int no_threads;

void central_barrier(bool& local_sense, size_t nr_threads)
{
    local_sense = !local_sense; // each processor toggles its own sense
    if (atomic_fetch_sub(&cnt, 1) == 1) {
        cnt = nr_threads;   // atomic store?
        sense.store(local_sense); // last processor toggles global sense 
    } else {
        while (local_sense != sense.load()) {}
    }
}

ThreadPool::~ThreadPool() 
{
    work_done_.store(true);
#ifdef USE_BOOST_BARRIER
    start_barrier_->wait();
    end_barrier_->wait();
#else
    central_barrier(local_sense, size_ + 1);
    central_barrier(local_sense, size_ + 1);
#endif

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
    if (size_)
        workers_.resize(size_);
    cnt = nr_threads + 1;
    local_sense = true;
    no_threads = nr_threads + 1;

#ifdef USE_BOOST_BARRIER
    start_barrier_ = make_shared<boost::barrier>(size_ + 1);
    end_barrier_ = make_shared<boost::barrier>(size_ + 1);
#endif

    for (size_t i = 0; i < size_; i++) {
        workers_[i].tid_ = i;
        workers_[i].local_sense = true;
#ifdef USE_BOOST_BARRIER
        workers_[i].start_barrier_ = start_barrier_;
        workers_[i].end_barrier_ = end_barrier_;
#endif
        // You can pass a member function pointer as the function, provided you
        // supply a suitable object pointer as the first argument
        workers_[i].thread_ = make_shared<boost::thread>
            (&ThreadPool::Run, this, workers_.data() + i);  // exception throwing try-catch
    }

#ifdef USE_BOOST_BARRIER
    start_barrier_->wait();
    end_barrier_->wait();
#else
    central_barrier(local_sense, nr_threads + 1);
#endif
}

void ThreadPool::Run(Worker *worker)
{
    tsc_t thread_tsc;
    tsc_init(&thread_tsc);

#ifdef USE_BOOST_BARRIER
    worker->start_barrier_->wait();
    worker->end_barrier_->wait();
#else
    central_barrier(worker->local_sense, no_threads);
#endif

    // Wait for the main thread to set a kernel to be executed
#ifdef USE_BOOST_BARRIER
    worker->start_barrier_->wait();
#else
    central_barrier(worker->local_sense, no_threads);
#endif
    while (!work_done_.load()) {
        switch (worker->job_) {
        case SPMV:
            tsc_start(&thread_tsc);
            do_matvec_thread(worker->data_);
            tsc_pause(&thread_tsc);
            break;
        default:
            break;
        }
#ifdef USE_BOOST_BARRIER
        worker->end_barrier_->wait();
        // Wait for a new kernel to be set
        worker->start_barrier_->wait();
#else
        central_barrier(worker->local_sense, no_threads);
        // Wait for a new kernel to be set
        central_barrier(worker->local_sense, no_threads);
#endif
    }
#ifdef USE_BOOST_BARRIER
    worker->end_barrier_->wait();
#else
    central_barrier(worker->local_sense, no_threads);
#endif

	spm_mt_thread_t *spm_thread = (spm_mt_thread_t *) worker->data_;
    spm_thread->secs = tsc_getsecs(&thread_tsc);
    tsc_shut(&thread_tsc);
}

void ThreadPool::SetKernel(int kernel_id)
{
    for (size_t i = 0; i < size_; i++) {
        workers_[i].SetJob(kernel_id);
    }
}
