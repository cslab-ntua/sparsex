/*
 * ThreadPool.hpp --  Thread management in SparseX.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include "CsxMatvec.hpp"

#include <atomic>
#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <vector>
#include "tsc.h"

#define SPMV 42

void central_barrier(bool& local_sense, size_t nr_threads);
class ThreadPool;

class Worker
{
public:
    Worker() {}
    Worker(size_t tid) : tid_(tid), job_(-1), local_sense(true) {}

    void SetJob(int job_id)
    {
        job_ = job_id;
    }

private:
    size_t tid_;
    int job_;
    bool local_sense;
    std::shared_ptr<boost::thread> thread_;  // use unique_ptr when C++14 available (make_unique)
    void *data_;
#ifdef USE_BOOST_BARRIER
    std::shared_ptr<boost::barrier> start_barrier_, end_barrier_;
#endif

    friend class ThreadPool;
};

class ThreadPool
{
public:
    static ThreadPool &GetInstance()
    {
        static ThreadPool instance;
        return instance;
    }

    ~ThreadPool();
	void InitThreads(size_t nr_threads);
    void Run(Worker *worker);
    void SetKernel(int kernel_id);

    void SetWorkerData(int worker_id, void *data)
    {
        workers_[worker_id].data_ = data;
    }

    size_t GetSize()
    {
        return size_;
    }

public:
#ifdef USE_BOOST_BARRIER
    std::shared_ptr<boost::barrier> start_barrier_;
    std::shared_ptr<boost::barrier> end_barrier_;
#endif
    bool local_sense;
private:
    size_t size_;
    atomic_bool work_done_;
    vector<Worker> workers_;

    ThreadPool() : local_sense(true), size_(0), work_done_(false) {}
    ThreadPool(ThreadPool const&);
    ThreadPool& operator=(ThreadPool const&);
};

#endif  // THREAD_POOL_HPP
