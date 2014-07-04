/*
 * ThreadPool.hpp --  Thread management in SparseX.
 *
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_THREAD_POOL_HPP
#define SPARSEX_INTERNALS_THREAD_POOL_HPP

#include "sparsex/internals/CsxMatvec.hpp"
#include "sparsex/internals/Barrier.hpp"

#include <atomic>
#include <boost/thread/thread.hpp>
#include <vector>

using namespace std;

enum Op {
    IDLE = -1,
    SPMV_MULT,
    SPMV_MULT_SYM,
    SPMV_KERNEL,
    SPMV_KERNEL_SYM
};

extern atomic<int> barrier_cnt;

class ThreadPool;

class Worker
{
public:
    Worker() : sense_(1), job_(IDLE) {}

    void SetId(size_t id)
    {
        id_ = id;
    }

    void SetJob(int job_id)
    {
        job_ = job_id;
    }

    void SetData(void *data)
    {
        data_ = data;
    }

    int GetJob()
    {
        return job_;
    }

    int *GetSense()
    {
        return &sense_;
    }
    
private:
    size_t id_;
    int sense_;
    int job_;
public:
    // use unique_ptr when C++14 available (make_unique)
    std::shared_ptr<boost::thread> thread_;
    void *data_;    // FIXME
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

    void SetWorkerData(size_t worker_id, void *data)
    {
        workers_[worker_id].SetData(data);
    }

    int *GetSense()
    {
        return &sense_;
    }

    size_t GetSize()
    {
        return size_;
    }

private:
    int sense_;
    size_t size_;
    atomic_bool work_done_;
    vector<Worker> workers_;

    ThreadPool() : sense_(1), size_(0), work_done_(false) {}
    ThreadPool(ThreadPool const&);
    ThreadPool& operator=(ThreadPool const&);
};

#endif  // SPARSEX_INTERNALS_THREAD_POOL_HPP
