/*
 * runtime.cc -- Front-end utilities for runtime configuration
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

#include "runtime.h"
#include "SparseUtil.h"
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
    mnemonic_map_.insert(MnemonicMap::value_type(RtNrThreads,
                                                 "libcsx.rt.nr_threads"));
    mnemonic_map_.insert(MnemonicMap::value_type(RtCpuAffinity,
                                                 "libcsx.rt.cpu_affinity"));
    mnemonic_map_.insert(MnemonicMap::value_type(PreprocXform,
                                                 "libcsx.preproc.xform"));
    mnemonic_map_.insert(MnemonicMap::value_type(PreprocMethod,
                                                 "libcsx.preproc.sampling"));
    mnemonic_map_.insert(MnemonicMap::value_type(
                             PreprocNrSamples,
                             "libcsx.preproc.sampling.nr_samples"));
    mnemonic_map_.insert(MnemonicMap::value_type(
                             PreprocSamplingPortion,
                             "libcsx.preproc.sampling.portion"));
    mnemonic_map_.insert(MnemonicMap::value_type(MatrixSymmetric,
                                                 "libcsx.matrix.symmetric"));
    mnemonic_map_.insert(MnemonicMap::value_type(MatrixSplitBlocks,
                                                 "libcsx.matrix.split_blocks"));
    mnemonic_map_.insert(MnemonicMap::value_type(MatrixFullColind,
                                                 "libcsx.matrix.full_colind"));
    mnemonic_map_.insert(MnemonicMap::value_type(
                             MatrixMinUnitSize,
                             "libcsx.matrix.min_unit_size"));
    mnemonic_map_.insert(MnemonicMap::value_type(
                             MatrixMaxUnitSize,
                             "libcsx.matrix.max_unit_size"));
    mnemonic_map_.insert(MnemonicMap::value_type(
                             MatrixMinCoverage,
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
        cout << "Window size: Not set\n";
    }

    const char *samples_str = getenv("SAMPLES");
    if (samples_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "portion");
        SetProperty(RuntimeConfiguration::PreprocNrSamples,
                    string(samples_str));
    } else {
        cout << "Number of samples: Not set\n";
    }

    const char *sampling_portion_str = getenv("SAMPLING_PORTION");
    if (sampling_portion_str) {
        // Automatically enable sampling
        SetProperty(RuntimeConfiguration::PreprocMethod, "portion");
        SetProperty(RuntimeConfiguration::PreprocSamplingPortion,
                    string(sampling_portion_str));
    } else {
        cout << "sampling is disabled\n";
    }
    
    return *this;
}

///> Max deltas that an encoding type may have. This is restricted by the
///  number of bits used for encoding the different patterns in CSX.
///
///  @see ctl_ll.h
#define DELTAS_MAX  CTL_PATTERNS_MAX

// malloc wrapper
#define xmalloc(x)                         \
({                                         \
    void *ret_;                            \
    ret_ = malloc(x);                      \
    if (ret_ == NULL){                     \
        cerr << __FUNCTION__          \
                  << " " << __FILE__       \
                  << ":" << __LINE__       \
                  << ": malloc failed\n";  \
        exit(1);                           \
    }                                      \
    ret_;                                  \
})

vector<size_t> &ParseOptionMT(string str, vector<size_t> &affinity)
{
    vector<string> mt_split;

    boost::split(mt_split, str, boost::algorithm::is_any_of(","),
                 boost::algorithm::token_compress_on);
    
    for (size_t i = 0; i < mt_split.size(); ++i) {
        affinity.push_back(boost::lexical_cast<size_t,string>(mt_split[i]));
    }
    
    // Printing
    cout << "MT_CONF=";
    for (size_t i = 0; i < affinity.size(); ++i) {
        if (i != 0)
            cout << ",";
        cout << affinity[i];
    }
    cout << "\n";
    return affinity;
}

static void ParseOptionXform(string str, int **xform_buf, int ***deltas)
{
    EncodingSequence encseq(str);
    cout << "Encseq: " << encseq << "\n";
    *xform_buf = (int *) xmalloc(XFORM_MAX * sizeof(int));
    *deltas = (int **) xmalloc(XFORM_MAX * sizeof(int *));
    for (size_t i = 0; i < XFORM_MAX; ++i)
        (*deltas)[i] = (int *) xmalloc(DELTAS_MAX * sizeof(int));

    vector<string> split_in_pairs;
    boost::algorithm::find_all_regex(split_in_pairs, str, 
                                     boost::regex
                                     ("[a-z]+(\\{[0-9]+(,[0-9]+)*\\})?"));
    vector<string> pair;
    vector<string> delta_tokens;
    size_t index = 0;
    size_t delta_index = 0;

    BOOST_FOREACH(string &token, split_in_pairs) {
        boost::split(pair, token, boost::algorithm::is_any_of("{"),
                     boost::algorithm::token_compress_on);

        // Find xform
        int option = XformIndex[pair[0]];
        if (option < XformBlockRow) {
            (*xform_buf)[index++] = option;
        } else if (option == XformBlockRow) {
            for (size_t j = BLOCK_TYPE_START; j < BLOCK_COL_START; ++j)
                (*xform_buf)[index++] = j;
        } else if (option == XformBlockCol) {
            for (size_t j = BLOCK_COL_START; j < BLOCK_TYPE_END; ++j)
                (*xform_buf)[index++] = j;
        } else {    // All
            for (size_t j = HORIZONTAL; j < XFORM_MAX; ++j)
                (*xform_buf)[index++] = j;
        }
        
        // Find deltas
        if (pair.size() > 1) {
            boost::algorithm::find_all_regex(delta_tokens, pair[1], 
                                             boost::regex("[0-9]+"));
            for (size_t j = 0; j < delta_tokens.size(); ++j) {
                (*deltas)[delta_index][j] = boost::lexical_cast<int,string>
                    (delta_tokens[j]);
            }
            (*deltas)[delta_index][delta_tokens.size()] = -1;
            delta_index++;
            if (option == XformBlockRow || option == XformBlockCol) {
                for (size_t k = 0; k < 7; ++k) {
                    (*deltas)[delta_index + k] = (*deltas)[delta_index - 1];
                }
                delta_index += 7;
            } else if (option == XformBlockAll) {
                for (size_t k = 0; k < 17; ++k) {
                    (*deltas)[delta_index + k] = (*deltas)[delta_index - 1];
                }
            }
        }
    }
    (*xform_buf)[index] = -1;

    // Printing
    cout << "Encoding type: ";
    for (size_t i = 0; (*xform_buf)[i] != -1; ++i) {
        if (i != 0)
            cout << ", ";
        cout << SpmTypesNames[(*xform_buf)[i]];
    }
    cout << "\n";   

    // In case there are no deltas free memory
    if (pair.size() == 1) {
        for (size_t i = 0; i < XFORM_MAX; ++i) {
            free((*deltas)[i]);
        }
        free(*deltas);
        (*deltas) = NULL;
    } else {
        cout << "Deltas to Encode: ";
        for (size_t i = 0; (*xform_buf)[i] != -1; ++i) {
            if (i != 0)
                cout << "}, ";
            cout << "{";
            assert((*deltas)[i][0] != -1);
            cout << (*deltas)[i][0];
            for (size_t j = 1; (*deltas)[i][j] != -1; ++j)
                cout << "," << (*deltas)[i][j];
        }
        cout << "}\n";
    }
}

void RuntimeContext::SetRuntimeContext(const RuntimeConfiguration &conf)
{
    cpu_affinity_ =
        ParseOptionMT(conf.GetProperty<string>(
                          RuntimeConfiguration::RtCpuAffinity),
                      cpu_affinity_);
}
