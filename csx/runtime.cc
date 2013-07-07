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
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

using namespace std;

#define BLOCK_ROW 5
#define BLOCK_COL 6
#define ALL       7
namespace csx {

static std::map<string, int> XformOpt = boost::assign::map_list_of
    ("none", 0)
    ("h", 1)
    ("v", 2)
    ("d", 3)
    ("ad", 4)
    ("br", 5)
    ("bc", 6)
    ("all", 7);

/**
 *  Load default CSX runtime properties
 */
RuntimeConfiguration::PropertyMap RuntimeConfiguration::DefaultProperties()
{
    return boost::assign::map_list_of
        (RtNrThreads, "1")
        (RtCpuAffinity, "0")
        (PreprocXform, "none")
        (PreprocSampling, "true")
        (PreprocNrSamples, "48")
        (PreprocSamplingPortion, "0.01")
        (PreprocWindowSize, "auto")
        (MatrixSymmetric, "false")
        (MatrixSplitBlocks, "true")
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
    mnemonic_map_.insert(MnemonicMap::value_type(PreprocSampling,
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
    if (mt_conf_str)
        SetProperty(RuntimeConfiguration::RtCpuAffinity, string(mt_conf_str));

    const char *xform_conf_str = getenv("XFORM_CONF");
    if (xform_conf_str)
        SetProperty(RuntimeConfiguration::PreprocXform, string(xform_conf_str));

    const char *wsize_str = getenv("WINDOW_SIZE");
    if (wsize_str) {
        SetProperty(RuntimeConfiguration::PreprocWindowSize, string(wsize_str));
    } else {
        cout << "Window size: Not set" << endl;
    }

    const char *samples_str = getenv("SAMPLES");
    if (samples_str) {
        SetProperty(RuntimeConfiguration::PreprocNrSamples,
                    string(samples_str));
    } else {
        cout << "Number of samples: Not set" << endl;
    }

    const char *sampling_portion_str = getenv("SAMPLING_PORTION");
    if (sampling_portion_str) {
        SetProperty(RuntimeConfiguration::PreprocSamplingPortion,
                    string(sampling_portion_str));
    } else {
        cout << "Sampling portion: Not set" << endl;
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

static void ParseOptionMT(string str, size_t **affinity, size_t &nr_threads)
{
    vector<string> mt_split;

    boost::split(mt_split, str, boost::algorithm::is_any_of(","),
                 boost::algorithm::token_compress_on);
    
    nr_threads = mt_split.size();
    *affinity = (size_t *) xmalloc(nr_threads * sizeof(size_t));
    
    for (size_t i = 0; i < nr_threads; i++) {
        (*affinity)[i] = boost::lexical_cast<size_t,string>(mt_split[i]);
    }
    
    // Printing
    cout << "MT_CONF=";
    for (size_t i = 0; i < nr_threads; i++) {
        if (i != 0)
            cout << ",";
        cout << (*affinity)[i];
    }
    cout << endl;
}

static void ParseOptionXform(string str, int **xform_buf, int ***deltas)
{
    *xform_buf = (int *) xmalloc(XFORM_MAX * sizeof(int));
    *deltas = (int **) xmalloc(XFORM_MAX * sizeof(int *));
    for (size_t i = 0; i < XFORM_MAX; i++) 
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
        int option = XformOpt[pair[0]];
        if (option < BLOCK_ROW)
            (*xform_buf)[index++] = option;
        else if (option == BLOCK_ROW) {
            for (size_t j = 7; j < 14; j++)
                (*xform_buf)[index++] = j;
        }
        else if (option == BLOCK_COL) {
            for (size_t j = 16; j < 23; j++)
                (*xform_buf)[index++] = j;
        } else {    // All
            for (size_t j = 1; j <= 4; j++)
                (*xform_buf)[index++] = j;
            for (size_t j = 7; j <= 13; j++)
                (*xform_buf)[index++] = j;
            for (size_t j = 16; j <= 22; j++)
                (*xform_buf)[index++] = j;
        }
        
        // Find deltas
        if (pair.size() > 1) {
            boost::algorithm::find_all_regex(delta_tokens, pair[1], 
                                             boost::regex("[0-9]+"));
            for (size_t j = 0; j < delta_tokens.size(); j++) {
                (*deltas)[delta_index][j] = boost::lexical_cast<int,string>
                    (delta_tokens[j]);
            }
            (*deltas)[delta_index][delta_tokens.size()] = -1;
            delta_index++;
            if (option == BLOCK_ROW || option == BLOCK_COL) {
                for (size_t k = 0; k < 7; k++) {
                    (*deltas)[delta_index + k] = (*deltas)[delta_index - 1];
                }
                delta_index += 7;
            } else if (option == ALL) {
                for (size_t k = 0; k < 17; k++) {
                    (*deltas)[delta_index + k] = (*deltas)[delta_index - 1];
                }
            }
        }
    }
    (*xform_buf)[index] = -1;

    // Printing
    cout << "Encoding type: ";
    for (size_t i = 0; (*xform_buf)[i] != -1; i++) {
        if (i != 0)
            cout << ", ";
        cout << SpmTypesNames[(*xform_buf)[i]];
    }
    cout << endl;   

    // In case there are no deltas free memory
    if (pair.size() == 1) {
        for (size_t i = 0; i < XFORM_MAX; i++) {
            free((*deltas)[i]);
        }
        free(*deltas);
        (*deltas) = NULL;
    } else {
        cout << "Deltas to Encode: ";
        for (size_t i = 0; (*xform_buf)[i] != -1; i++) {
            if (i != 0)
                cout << "}, ";
            cout << "{";
            assert((*deltas)[i][0] != -1);
            cout << (*deltas)[i][0];
            for (size_t j = 1; (*deltas)[i][j] != -1; j++)
                cout << "," << (*deltas)[i][j];
        }
        cout << "}" << endl;
    }
}

static void ParseOptionWindowSize(string str, size_t &wsize)
{
    wsize = boost::lexical_cast<size_t, string>(str);
    cout << "Window size: " << wsize << endl;
}

static void ParseOptionSamples(string str, size_t &samples_max)
{
    if (str == "__MAX__")
        samples_max = std::numeric_limits<uint64_t>::max();
    else
        samples_max = boost::lexical_cast<size_t, string>(str);
    cout << "Number of samples: " << samples_max << endl;
}

static void ParseOptionPortion(string str, double &sampling_portion)
{
    sampling_portion = boost::lexical_cast<double, string>(str);
    cout << "Sampling portion: " << sampling_portion << endl;
}

static void ParseOptionSymmetric(string str, bool &symmetric)
{
    (str == "true") ? symmetric = true : symmetric = false;
}

static void ParseOptionSplitBlocks(string str, bool &split_blocks)
{
    (str == "true") ? split_blocks = true : split_blocks = false;   
} 

}

void RuntimeContext::SetRuntimeContext(const RuntimeConfiguration &conf)
{
    ParseOptionMT(conf.GetProperty(RuntimeConfiguration::RtCpuAffinity),
                  &affinity_, nr_threads_);

    // Initialize the CSX JIT execution engine
    engine_ = &CsxJitInit();
}

void CsxContext::SetCsxContext(const RuntimeConfiguration &conf)
{
    ParseOptionSymmetric(conf.GetProperty(
                             RuntimeConfiguration::MatrixSymmetric),
                         symmetric_);
    ParseOptionSplitBlocks(conf.GetProperty(
                               RuntimeConfiguration::MatrixSplitBlocks),
                           split_blocks_);
    ParseOptionXform(conf.GetProperty(RuntimeConfiguration::PreprocXform),
                     &xform_buf_, &deltas_);
    ParseOptionWindowSize(conf.GetProperty(
                              RuntimeConfiguration::PreprocWindowSize),
                          wsize_); 
    ParseOptionSamples(conf.GetProperty(RuntimeConfiguration::PreprocNrSamples),
                       samples_max_);
    ParseOptionPortion(conf.GetProperty(
                           RuntimeConfiguration::PreprocSamplingPortion),
                       sampling_portion_);
}
