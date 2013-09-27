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

#include "csx.h"
#include "jit.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sched.h>
#include <string>
#include <boost/bimap.hpp>
#include <boost/lexical_cast.hpp>

extern "C" {
#include "spm_mt.h"
}

using namespace std;

namespace csx {

template<typename IndexType, typename ValueType>
class SparseInternal;

/* Runtime configuration class */
class RuntimeConfiguration
{
public:
    enum Property {
        RtNrThreads,
        RtCpuAffinity,
        PreprocXform,
        PreprocMethod,
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
        LogFilename,
        LogLevel,
    };

    static RuntimeConfiguration &GetInstance()
    {
        static RuntimeConfiguration instance;
        return instance;
    }

    RuntimeConfiguration &LoadFromEnv();

    template<typename Target>
    Target GetProperty(const Property &key) const
    {
        PropertyMap::const_iterator iter = property_map_.find(key);
        return boost::lexical_cast<Target, string>(iter->second);
    }

    void SetProperty(const Property &key, const string &value)
    {
        property_map_[key] = value;
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

// Special conversions
// Support `true' and `false' in boolean conversions
template<>
bool RuntimeConfiguration::GetProperty(const Property &key) const
{
    PropertyMap::const_iterator iter = property_map_.find(key);
    assert(iter != property_map_.end());
    istringstream ss(iter->second);
    bool ret;
    ss >> std::boolalpha >> ret;
    return ret;
}

template<>
PreprocessingMethod RuntimeConfiguration::GetProperty(const Property &key) const
{
    PropertyMap::const_iterator iter = property_map_.find(key);
    PreprocessingMethod ret(iter->second);
    return ret;
}

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
        return cpu_affinity_.size();
    }

    size_t GetAffinity(size_t tid) const
    {
        return cpu_affinity_[tid];
    }

    CsxExecutionEngine &GetEngine() const
    {
        return *engine_;
    }

private:
    RuntimeContext()
        : engine_(0)
    {
        // Initialize the CSX JIT execution engine
        engine_ = &CsxJitInit();
    }

    ~RuntimeContext()
    {
        // do not free engine_; we do not own it
    }

    RuntimeContext(const RuntimeContext &); // Do not implement
    RuntimeContext &operator=(const RuntimeContext &); // Do not implement

    vector<size_t> cpu_affinity_;
    CsxExecutionEngine *engine_;
};

/* Class responsible for holding per thread information */
template<class InternalType>
class ThreadContext
{
public:
    typedef typename InternalType::idx_t idx_t;
    typedef typename InternalType::val_t val_t;

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

    void SetData(SparseInternal<idx_t, val_t> *spms, spm_mt_t *spm_mt);
    void SetDataSym(SparsePartitionSym<idx_t, val_t> *spms,
                    spm_mt_t *spm_mt);
 
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

    CsxManager<idx_t, val_t> *GetCsxManager()
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
    CsxManager<idx_t, val_t> *csxmg_;
    ostringstream *buffer_;
    double sampling_prob_;
};

template<class InternalType>
void ThreadContext<InternalType>::
SetData(SparseInternal<idx_t, val_t> *spms, spm_mt_t *spm_mt)
{
    spm_ = spms->GetPartition(id_);
    csxmg_ = new CsxManager<idx_t, val_t>(spm_);
    spm_encoded_ = &spm_mt->spm_threads[id_]; 
    buffer_ = new ostringstream("");
#ifdef SPM_NUMA
    // Enable the full-column-index optimization for NUMA architectures
    csxmg_->SetFullColumnIndices(true);
#endif
}

template<class InternalType>
void ThreadContext<InternalType>::
SetDataSym(SparsePartitionSym<idx_t, val_t> *spms, spm_mt_t *spm_mt)
{
    spm_ = spms + id_;
    csxmg_ = new CsxManager<idx_t, val_t>(spms + id_);
    spm_encoded_ = &spm_mt->spm_threads[id_]; 
    buffer_ = new ostringstream("");
#ifdef SPM_NUMA
    // Enable the full-column-index optimization for NUMA architectures
    csxmg_->SetFullColumnIndices(true);
#endif
}

}   // end of namespace csx

#endif // RUNTIME_H__
