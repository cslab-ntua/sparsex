/* -*- C++ -*-
 * 
 * runtime.h -- Front-end utilities for runtime configuration
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2012-2013, Athena Elafrou
 * Copyright (C) 2013,      Vasileios Karakasis
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
#include <boost/bimap.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>

#include "SparseInternal.h"
#include "csx.h"
#include "jit.h"

extern "C" {
#include "mt_lib.h" //To be removed
#include "spm_mt.h"
}

using namespace std;

namespace csx {

/* Runtime configuration class */
class RuntimeConfiguration
{
public:
    enum Property {
        RtNrThreads,
        RtCpuAffinity,
        PreprocXform,
        PreprocSampling,
        PreprocNrSamples,
        PreprocSamplingPortion,
        PreprocWindowSize,
        MatrixSymmetric,
        MatrixSplitBlocks,
        MatrixFullColind,
        MatrixOneDimBlocks,
        MatrixMinUnitSize,
        MatrixMaxUnitSize,
        MatrixMinCoverage,
    };

    static RuntimeConfiguration &GetInstance()
    {
        static RuntimeConfiguration instance;
        return instance;
    }

    RuntimeConfiguration &LoadFromEnv();

    const string &GetProperty(const Property &key) const
    {
        PropertyMap::const_iterator iter = property_map_.find(key);
        return iter->second;
    }

    void SetProperty(const Property &key, const string &value)
    {
        property_map_.insert(make_pair(key, value));
    }

private:
    typedef std::map<Property, string> PropertyMap;
    typedef boost::bimap<Property, string> MnemonicMap;

    RuntimeConfiguration()
        : property_map_(DefaultProperties())
    {
        InitPropertyMnemonics();
    }
    
    ~RuntimeConfiguration() {}

    /**
     *  Load default CSX runtime properties
     */
    PropertyMap DefaultProperties();
    void InitPropertyMnemonics();

    PropertyMap property_map_;
    MnemonicMap mnemonic_map_;
};

/* Singleton class holding runtime information */
class RuntimeContext
{
public:

    static RuntimeContext &GetInstance()
    {
        static RuntimeContext instance;
        return instance;
    }

    void SetRuntimeContext(const RuntimeConfiguration &conf);

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
    RuntimeContext()
        : nr_threads_(1), 
          affinity_(NULL),
          engine_(NULL)
    {}

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
    static CsxContext &GetInstance()
    {
        static CsxContext instance;
        return instance;
    }

    void SetCsxContext(const RuntimeConfiguration &conf);

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
    CsxContext()
        : symmetric_(false),
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

    bool symmetric_, split_blocks_;
    int *xform_buf_;
    int **deltas_;
    size_t wsize_;
    size_t samples_max_;
    double sampling_portion_;
};

/* Class responsible for holding per thread information */
template<class InternalType>
class ThreadContext
{
public:
    typedef typename InternalType::index_t index_t;
    typedef typename InternalType::value_t value_t;

    ThreadContext() : 
        id_(0), 
        cpu_(0), 
        node_(0), 
        spm_(NULL), 
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

    //template<typename IndexType, typename ValueType>
    void SetData(csx::SparseInternal<index_t, value_t> *spms, spm_mt_t *spm_mt);

#if 0   // SYM
    void SetDataSym(SparsePartitionSym<index_t, value_t> *spms,
                    spm_mt_t *spm_mt);
#endif  // SYM    
 
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

    InternalType *GetSpm()
    {
        return spm_;
    }

    spm_mt_thread_t *GetSpmEncoded()
    {
        return spm_encoded_;
    }

    CsxManager<index_t, value_t> *GetCsxManager()
    {
        return csxmg_;
    }

    ostringstream& GetBuffer()
    {
        return *buffer_;
    }

private:
    size_t id_, cpu_, node_;
    InternalType *spm_;
    spm_mt_thread_t *spm_encoded_;
    CsxManager<index_t, value_t> *csxmg_;
    ostringstream *buffer_;
    double sampling_prob_;
};

template<class InternalType>
void ThreadContext<InternalType>::
SetData(SparseInternal<index_t, value_t> *spms, spm_mt_t *spm_mt)
{
    spm_ = spms->GetPartition(id_);
    csxmg_ = new CsxManager<index_t, value_t>(spm_);
    spm_encoded_ = &spm_mt->spm_threads[id_]; 
    buffer_ = new ostringstream("");
#ifdef SPM_NUMA
    // Enable the full-column-index optimization for NUMA architectures
    csxmg_->SetFullColumnIndices(true);
#endif
}

#if 0   // SYM

template<class InternalType>
void ThreadContext<InternalType>::
SetDataSym(SparsePartitionSym<index_t, value_t> *spms, spm_mt_t *spm_mt)
{
    spm_ = spms + id_;
    csxmg_ = new CsxManager<index_t, value_t>(spms + id_);
    spm_encoded_ = &spm_mt->spm_threads[id_]; 
    buffer_ = new ostringstream("");
#ifdef SPM_NUMA
    // Enable the full-column-index optimization for NUMA architectures
    csxmg_->SetFullColumnIndices(true);
#endif
}

#endif  // SYM

}   // end of namespace csx

#endif // RUNTIME_H__
