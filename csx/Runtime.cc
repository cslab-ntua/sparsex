/*
 * Runtime.cc -- Front-end utilities for runtime configuration
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

#include "Runtime.hpp"
#include "SparseUtil.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/assign/list_of.hpp>

using namespace std;
using namespace boost;

static std::map<string, int> XformIndex = boost::assign::map_list_of
    ("none", 0)
    ("h", 1)
    ("v", 2)
    ("d", 3)
    ("ad", 4)
    ("br", 5)
    ("bc", 6)
    ("all", 7);

const int XformBlockRow = XformIndex["br"];
const int XformBlockCol = XformIndex["bc"];
const int XformBlockAll = XformIndex["all"];

/**
 *  Load default CSX runtime properties
 */
RuntimeConfiguration::PropertyMap RuntimeConfiguration::DefaultProperties()
{
    return boost::assign::map_list_of
        (RtNrThreads, "1")
        (RtCpuAffinity, "0")
        (PreprocXform, "none")
        (PreprocMethod, "none")
        (PreprocNrSamples, "10")
        (PreprocSamplingPortion, "0.01")
        (PreprocWindowSize, "0")    // change to sth more meaningful
        (MatrixSymmetric, "false")
        (MatrixSplitBlocks, "true")
        (MatrixOneDimBlocks, "false")
#ifdef SPM_NUMA
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
    mnemonic_map_.insert(mnemonic_type(RtNrThreads, 
                                       "libcsx.rt.nr_threads"));
    mnemonic_map_.insert(mnemonic_type(RtCpuAffinity, 
                                       "libcsx.rt.cpu_affinity"));
    mnemonic_map_.insert(mnemonic_type(PreprocXform, 
                                       "libcsx.preproc.xform"));
    mnemonic_map_.insert(mnemonic_type(PreprocMethod, 
                                       "libcsx.preproc.sampling"));
    mnemonic_map_.insert(mnemonic_type(PreprocNrSamples,
                                       "libcsx.preproc.sampling.nr_samples"));
    mnemonic_map_.insert(mnemonic_type(PreprocSamplingPortion,
                                       "libcsx.preproc.sampling.portion"));
    mnemonic_map_.insert(mnemonic_type(MatrixSymmetric,
                                       "libcsx.matrix.symmetric"));
    mnemonic_map_.insert(mnemonic_type(MatrixSplitBlocks, 
                                       "libcsx.matrix.split_blocks"));
    mnemonic_map_.insert(mnemonic_type(MatrixFullColind,
                                       "libcsx.matrix.full_colind"));
    mnemonic_map_.insert(mnemonic_type(MatrixMinUnitSize,
                                       "libcsx.matrix.min_unit_size"));
    mnemonic_map_.insert(mnemonic_type(MatrixMaxUnitSize,
                                       "libcsx.matrix.max_unit_size"));
    mnemonic_map_.insert(mnemonic_type(MatrixMinCoverage,
                                       "libcsx.matrix.min_coverage"));
}

RuntimeConfiguration &RuntimeConfiguration::LoadFromEnv()
{
    const char *mt_conf_str = getenv("MT_CONF");
    if (mt_conf_str) {
        SetProperty(RuntimeConfiguration::RtCpuAffinity, string(mt_conf_str));
    }

    const char *xform_conf_str = getenv("XFORM_CONF");
    if (xform_conf_str)
        SetProperty(RuntimeConfiguration::PreprocXform, string(xform_conf_str));

    const char *wsize_str = getenv("WINDOW_SIZE");
    if (wsize_str) {
        SetProperty(RuntimeConfiguration::PreprocWindowSize, string(wsize_str));
    } else {
        // cout << "Window size: Not set\n";
    }

    const char *samples_str = getenv("SAMPLES");
    if (samples_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "portion");
        SetProperty(RuntimeConfiguration::PreprocNrSamples,
                    string(samples_str));
    } else {
        // cout << "Number of samples: Not set\n";
    }

    const char *sampling_portion_str = getenv("SAMPLING_PORTION");
    if (sampling_portion_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "portion");
        SetProperty(RuntimeConfiguration::PreprocSamplingPortion,
                    string(sampling_portion_str));
    } else {
        // cout << "sampling is disabled\n";
    }
    
    return *this;
}

vector<size_t> &ParseOptionMT(string str, vector<size_t> &affinity)
{
    vector<string> mt_split;

    boost::split(mt_split, str, boost::algorithm::is_any_of(","),
                 boost::algorithm::token_compress_on);
    
    for (size_t i = 0; i < mt_split.size(); ++i) {
        affinity.push_back(boost::lexical_cast<size_t,string>(mt_split[i]));
    }
    
    // Printing
    // cout << "MT_CONF=";
    // for (size_t i = 0; i < affinity.size(); ++i) {
    //     if (i != 0)
    //         cout << ",";
    //     cout << affinity[i];
    // }
    // cout << "\n";

    return affinity;
}

void RuntimeContext::SetRuntimeContext(const RuntimeConfiguration &conf)
{
    cpu_affinity_ =
        ParseOptionMT(conf.GetProperty<string>(
                          RuntimeConfiguration::RtCpuAffinity),
                      cpu_affinity_);
}
