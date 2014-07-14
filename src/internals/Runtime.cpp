/*
 * \file Runtime.cpp
 *
 * \brief Front-end utilities for runtime configuration
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2013, Vasileios Karakasis
 * Copyright (C) 2012-2013, Athena Elafrou
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Runtime.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/assign/list_of.hpp>

using namespace std;
using namespace boost;

namespace sparsex {
namespace runtime {

/**
 *  Load default CSX runtime properties
 */
RuntimeConfiguration::PropertyMap RuntimeConfiguration::DefaultProperties()
{
    return boost::assign::map_list_of
        (RtNrThreads, "1")
        (RtCpuAffinity, "0")
#if SPX_USE_NUMA
        (PreprocHeuristic, "cost")
#else
        (PreprocHeuristic, "ratio")
#endif
        (PreprocXform, "none")
        (PreprocMethod, "none")
        (PreprocNrSamples, "10")
        (PreprocSamplingPortion, "0.01")
        (PreprocWindowSize, "0")
        (MatrixSymmetric, "false")
        (MatrixSplitBlocks, "true")
        (MatrixOneDimBlocks, "false")
#if SPX_USE_NUMA
        (MatrixFullColind, "true")
#else
        (MatrixFullColind, "false")
#endif            
        (MatrixMinUnitSize, "4")
        (MatrixMaxUnitSize, "255")
        (MatrixMinCoverage, "0.1");
}

void RuntimeConfiguration::InitPropertyMnemonics()
{
    mnemonic_map_.insert(mnemonic_pair(RtNrThreads, 
                                       "spx.rt.nr_threads"));
    mnemonic_map_.insert(mnemonic_pair(RtCpuAffinity, 
                                       "spx.rt.cpu_affinity"));
    mnemonic_map_.insert(mnemonic_pair(PreprocHeuristic, 
                                       "spx.preproc.heuristic"));
    mnemonic_map_.insert(mnemonic_pair(PreprocXform, 
                                       "spx.preproc.xform"));
    mnemonic_map_.insert(mnemonic_pair(PreprocMethod, 
                                       "spx.preproc.sampling"));
    mnemonic_map_.insert(mnemonic_pair(PreprocNrSamples,
                                       "spx.preproc.sampling.nr_samples"));
    mnemonic_map_.insert(mnemonic_pair(PreprocSamplingPortion,
                                       "spx.preproc.sampling.portion"));
    mnemonic_map_.insert(mnemonic_pair(PreprocWindowSize,
                                       "spx.preproc.sampling.window_size"));
    mnemonic_map_.insert(mnemonic_pair(MatrixSymmetric,
                                       "spx.matrix.symmetric"));
    mnemonic_map_.insert(mnemonic_pair(MatrixSplitBlocks, 
                                       "spx.matrix.split_blocks"));
    mnemonic_map_.insert(mnemonic_pair(MatrixFullColind,
                                       "spx.matrix.full_colind"));
    mnemonic_map_.insert(mnemonic_pair(MatrixMinUnitSize,
                                       "spx.matrix.min_unit_size"));
    mnemonic_map_.insert(mnemonic_pair(MatrixMaxUnitSize,
                                       "spx.matrix.max_unit_size"));
    mnemonic_map_.insert(mnemonic_pair(MatrixMinCoverage,
                                       "spx.matrix.min_coverage"));
}

RuntimeConfiguration &RuntimeConfiguration::LoadFromEnv()
{
    const char *mt_conf_str = getenv("MT_CONF");
    if (mt_conf_str) {
        SetProperty(RuntimeConfiguration::RtCpuAffinity, string(mt_conf_str));
    }

    const char *xform_conf_str = getenv("XFORM_CONF");
    if (xform_conf_str) {
        SetProperty(RuntimeConfiguration::PreprocXform, string(xform_conf_str));
    }

    const char *wsize_str = getenv("WINDOW_SIZE");
    if (wsize_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "window");
        SetProperty(RuntimeConfiguration::PreprocWindowSize, string(wsize_str));
    }

    const char *samples_str = getenv("SAMPLES");
    if (samples_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "portion");
        SetProperty(RuntimeConfiguration::PreprocNrSamples,
                    string(samples_str));
    }

    const char *sampling_portion_str = getenv("SAMPLING_PORTION");
    if (sampling_portion_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "portion");
        SetProperty(RuntimeConfiguration::PreprocSamplingPortion,
                    string(sampling_portion_str));
    }

    return *this;
}

vector<size_t> &ParseOptionMT(string str, vector<size_t> &affinity)
{
    vector<string> mt_split;
    boost::split(mt_split, str, boost::algorithm::is_any_of(","),
                 boost::algorithm::token_compress_on);

    // Clear in case it has been already set.
    affinity.clear();
    for (size_t i = 0; i < mt_split.size(); ++i) {
        affinity.push_back(boost::lexical_cast<size_t, string>(mt_split[i]));
    }
    
    return affinity;
}

void RuntimeContext::CheckParams(const RuntimeConfiguration &conf)
{
    vector<size_t> affinity;
    affinity = ParseOptionMT(
        conf.GetProperty<string>(RuntimeConfiguration::RtCpuAffinity), affinity);

    LOG_INFO << "Number of threads: " << affinity.size() << "\n";
    stringstream os;
    os << "Thread affinity: {";
    for (size_t i = 0; i < affinity.size(); ++i) {
        if (i != 0)
            os << ",";
        os << affinity[i];
    }
    os << "}\n";
    LOG_INFO << os.str();
}

void RuntimeContext::SetRuntimeContext(const RuntimeConfiguration &conf)
{
    cpu_affinity_ =
        ParseOptionMT(conf.GetProperty<string>(
                          RuntimeConfiguration::RtCpuAffinity),
                      cpu_affinity_);
}

} // end of namespace runtime
} // end of namespace sparsex
