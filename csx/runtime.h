/* -*- C++ -*-
 * 
 * runtime.h -- Front-end utilities for runtime configuration
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef RUNTIME_H__
#define RUNTIME_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sched.h>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>

#include "spm.h"
#include "csx.h"
#include "jit.h"

extern "C" {
#include "mt_lib.h" //To be removed
#include "spm_mt.h"
}

using namespace std;
using namespace csx;

/* Runtime configuration class */
class Configuration
{
public:
    Configuration() {}
    ~Configuration() {}

    template<typename Field>
    string GetProperty(Field field) const
    {
        return get(conf_, field);
    }

    template<typename Field>
    void SetProperty(Field field, string new_value)
    {
        put(conf_, field, new_value);
    }

    enum ConfigInfo {
        /* User-defined */
        NrThreads,
        Affinity,
        Xform,
        WindowSize,
        Samples,
        SamplingPortion,
        Symmetric,
        SplitBlocks, 
        /* Implementation-defined */
        OneDimBlocks,
        FullColumnIndices,
        MinLimit,
        MaxLimit,
        MinPercentage
    };

private:
    typedef boost::unordered_map<ConfigInfo, string> Map;
    typedef boost::associative_property_map<Map> PropertyMap;
    static Map map_;
    static PropertyMap conf_;
};

/* 
 * Loads runtime configuration set by the user
 * in the global property map "config_map"
 */
Configuration &ConfigFromEnv(Configuration &conf);

/* Singleton class holding runtime information */
class RuntimeContext
{
public:

    static RuntimeContext &GetInstance()
    {
        static RuntimeContext instance;
        return instance;
    }

    void SetRuntimeContext(const Configuration &conf);

    size_t GetNrThreads() const
    {
        return nr_threads_;
    }

    size_t GetAffinity(size_t tid) const
    {
        return affinity_[tid];
    }

    CsxExecutionEngine &GetEngine() const
    {
        return *engine_;
    }

private:

    RuntimeContext() : 
        nr_threads_(1), 
        affinity_(NULL),
        engine_(NULL) {}

    ~RuntimeContext()
    {
        free(affinity_);
    }

    RuntimeContext(const RuntimeContext &); // Do not implement
    RuntimeContext &operator=(const RuntimeContext &); // Do not implement

    size_t nr_threads_, *affinity_;
    CsxExecutionEngine *engine_;    //smart pointer???
};

/* Class holding CSX-related configuration */
class CsxContext
{
public:

    CsxContext() :
        symmetric_(false),
        split_blocks_(true),
        xform_buf_(NULL), 
        deltas_(NULL),
        wsize_(0),
        samples_max_(0), 
        sampling_portion_(0) {}

    ~CsxContext()
    {
        if (xform_buf_)
            free(xform_buf_);
        if (deltas_) {
            for (size_t i = 0; i < XFORM_MAX; i++)
                free(deltas_[i]);
            free(deltas_);
        }
    }

    void SetCsxContext(const Configuration &conf);

    bool IsSymmetric() const
    {
        return symmetric_;
    }

    bool IsSplitBlocks() const
    {
        return split_blocks_;
    }

    int *GetXform() const   //add const
    {
        return xform_buf_;
    }
    
    int **GetDeltas() const
    {
        return deltas_;
    }

    size_t GetWindowSize() const
    {
        return wsize_;
    }

    size_t GetMaxSamples() const
    {
        return samples_max_;
    }

    double GetSamplingPortion() const
    {
        return sampling_portion_;
    }

private:
    bool symmetric_, split_blocks_;
    int *xform_buf_;
    int **deltas_;
    size_t wsize_;
    size_t samples_max_;
    double sampling_portion_;
};

/* Class responsible for holding per thread information */
class ThreadContext
{
public:

    ThreadContext() : 
        id_(0), 
        cpu_(0), 
        node_(0), 
        spm_(NULL), 
        spm_sym_(NULL),
        spm_encoded_(NULL),
        csxmg_(NULL),
        buffer_(NULL),
        sampling_prob_(0) {}

    ~ThreadContext() {}

    void SetId(size_t id)
    {
        id_ = id;
    }

    void SetCpu(size_t cpu)
    {
        cpu_ = cpu;
    }

    void SetNode(size_t node)
    {
        node_ = node;
    }

    void SetData(SPM* spms, SPMSym *spms_sym, spm_mt_t *spm_mt, bool symmetric);
 
    size_t GetId() const
    {
        return id_;
    }

    size_t GetCpu() const
    {
        return cpu_;
    }

    size_t GetNode() const
    {
        return node_;
    }

    SPM *GetSpm()
    {
        return spm_;
    }

    SPMSym *GetSpmSym()
    {
        return spm_sym_;
    }

    spm_mt_thread_t *GetSpmEncoded()
    {
        return spm_encoded_;
    }

    CsxManager *GetCsxManager()
    {
        return csxmg_;
    }

    ostringstream &GetBuffer()
    {
        return *buffer_;
    }

private:
    size_t id_, cpu_, node_;
    SPM *spm_;
    SPMSym *spm_sym_;
    spm_mt_thread_t *spm_encoded_;
    CsxManager *csxmg_; //smart pointer???
    ostringstream *buffer_;
    double sampling_prob_;
};

#endif // RUNTIME_H__
