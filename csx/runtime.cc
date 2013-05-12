/*
 * runtime.cc -- Front-end utilities for runtime configuration
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "runtime.h"
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

using namespace std;
using namespace csx;

#define BLOCK_ROW 5
#define BLOCK_COL 6
#define ALL       7

static std::map<string, int> XformOpt = boost::assign::map_list_of
    ("none", 0)
    ("h", 1)
    ("v", 2)
    ("d", 3)
    ("ad", 4)
    ("br", 5)
    ("bc", 6)
    ("all", 7);

// Initialization of local static data members
Configuration::Map Configuration::map_ = boost::assign::map_list_of
    (NrThreads, "1")
    (Affinity, "0")
    (Xform, "0")
    (WindowSize, "0")
    (Samples, "__MAX__")    
    (SamplingPortion, "0.0")
    (Symmetric, "false")
    (SplitBlocks, "true")
    (OneDimBlocks, "false")
#ifdef SPM_NUMA
    (FullColumnIndices, "true")
#else
    (FullColumnIndices, "false")
#endif
    (MinLimit, "4")
    (MaxLimit, "__MAX__")  //default -> long max
    (MinPercentage, "0.1");

Configuration::PropertyMap Configuration::conf_(Configuration::map_);

Configuration &ConfigFromEnv(Configuration &conf)
{
    const char *mt_conf_str = getenv("MT_CONF");
    if (mt_conf_str)
        conf.SetProperty(Configuration::Affinity, string(mt_conf_str));

    const char *xform_conf_str = getenv("XFORM_CONF");
    if (xform_conf_str)
        conf.SetProperty(Configuration::Xform, string(xform_conf_str));

    const char *wsize_str = getenv("WINDOW_SIZE");
    if (wsize_str) {
        conf.SetProperty(Configuration::WindowSize, string(wsize_str));
    } else {
        cout << "Window size: Not set" << endl;
    }

    const char *samples_str = getenv("SAMPLES");
    if (samples_str) {
        conf.SetProperty(Configuration::Samples, string(samples_str));
    } else {
        cout << "Number of samples: Not set" << endl;
    }

    const char *sampling_portion_str = getenv("SAMPLING_PORTION");
    if (sampling_portion_str) {
        conf.SetProperty(Configuration::SamplingPortion,
                         string(sampling_portion_str));
    } else {
        cout << "Sampling portion: Not set" << endl;
    }
    
    return conf;
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
                                     boost::regex("[a-z]+(\\{[0-9]+(,[0-9]+)*\\})?"));

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
                (*deltas)[delta_index][j] = boost::lexical_cast<int,string>(delta_tokens[j]);
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

    // In case there are no deltas release memory
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

void RuntimeContext::SetRuntimeContext(const Configuration &conf)
{
    ParseOptionMT(conf.GetProperty(Configuration::Affinity), &affinity_, nr_threads_);

    // Initialize the CSX JIT execution engine
    engine_ = &CsxJitInit();
}

void CsxContext::SetCsxContext(const Configuration &conf)
{
    ParseOptionSymmetric(conf.GetProperty(Configuration::Symmetric), symmetric_);
    ParseOptionSplitBlocks(conf.GetProperty(Configuration::SplitBlocks),
                           split_blocks_);
    ParseOptionXform(conf.GetProperty(Configuration::Xform), &xform_buf_, &deltas_);
    ParseOptionWindowSize(conf.GetProperty(Configuration::WindowSize), wsize_); 
    ParseOptionSamples(conf.GetProperty(Configuration::Samples), samples_max_);
    ParseOptionPortion(conf.GetProperty(Configuration::SamplingPortion),
                       sampling_portion_);
}

void ThreadContext::SetData(SPM* spms, SPMSym *spms_sym, spm_mt_t *spm_mt,
                            bool symmetric)
{
    if (!symmetric) {
        spm_ = spms + id_;
        spm_sym_ = NULL;
        csxmg_ = new CsxManager(spms + id_); //changed
    } else {
        spm_ = NULL;
        spm_sym_ = spms_sym + id_;
        csxmg_ = new CsxManager(spms_sym + id_);
    }
    spm_encoded_ = &spm_mt->spm_threads[id_]; 
    buffer_ = new ostringstream("");
#ifdef SPM_NUMA
    // Enable the full-column-index optimization for NUMA architectures
    csxmg_->SetFullColumnIndices(true);
#endif
}

/*static int SplitString(char *str, char **str_buf, const char *start_sep,
                       const char *end_sep)
{
    char *token = strtok(str, start_sep);
    int next = 0;
    int str_length = strcspn(token, end_sep);

    str_buf[next] = (char *) xmalloc((str_length + 1) * sizeof(char));
    strncpy(str_buf[next], token, str_length);
    str_buf[next][str_length] = 0;
    ++next;
    while ((token = strtok(NULL, start_sep)) != NULL) {
        str_length = strcspn(token, end_sep);
        str_buf[next] = (char *) xmalloc((str_length - 1) * sizeof(char));
        strncpy(str_buf[next], token, str_length);
        str_buf[next][str_length] = 0;
        ++next;
    }

    return next;
}

static void GetOptionEncodeDeltas(int ***deltas)
{
    char *encode_deltas_env = getenv("ENCODE_DELTAS");
    if (encode_deltas_env && strlen(encode_deltas_env)) {
        // Copy environment variable to avoid being altered from strtok()
        char *encode_deltas_str = (char *) xmalloc(strlen(encode_deltas_env)+1);
        strncpy(encode_deltas_str, encode_deltas_env,
                strlen(encode_deltas_env)+1);

        // Init matrix deltas.
        *deltas = (int **) xmalloc(XFORM_MAX * sizeof(int *));

        for (int i = 0; i < XFORM_MAX; i++) {
            (*deltas)[i] = (int *) xmalloc(DELTAS_MAX * sizeof(int));
        }

        // Fill deltas with the appropriate data.
        char **temp = (char **) xmalloc(XFORM_MAX * sizeof(char *));
        int temp_size = SplitString(encode_deltas_str, temp, "{", "}");

        for (int i = 0; i < temp_size; i++) {
            int j = 0;
            char *token = strtok(temp[i], ",");

            (*deltas)[i][j] = atoi(token);
            ++j;
            while ((token = strtok(NULL, ",")) != NULL) {
                (*deltas)[i][j] = atoi(token);
                ++j;
            }

            (*deltas)[i][j] = -1;
            free(temp[i]);
        }

        free(temp);

        // Print deltas
        cout << "Deltas to Encode: ";
        for (int i = 0; i < temp_size; i++) {
            if (i != 0)
                cout << "}, ";
            cout << "{";
            assert((*deltas)[i][0] != -1);
            cout << (*deltas)[i][0];
            for (int j = 1; (*deltas)[i][j] != -1; j++)
                cout << "," << (*deltas)[i][j];
        }

        cout << "}" << endl;
        free(encode_deltas_str);
    }
}

static void ParseOptionDeltas(string str, int ***deltas)
{
    if (!str.empty()) {
        vector<string> deltas_split;
        boost::algorithm::find_all_regex(deltas_split, str,
                                         boost::regex("[0-9]+(,[0-9]+)*"));

        // Init matrix deltas.
        *deltas = (int **) xmalloc((deltas_split.size()) * sizeof(int *));

        for (size_t i = 0; i < deltas_split.size(); i++) {
            (*deltas)[i] = (int *) xmalloc(DELTAS_MAX * sizeof(int));
        }

        for (size_t i = 0; i < deltas_split.size(); i++) {
            vector<string> tokens;
            boost::split(tokens, deltas_split[i], boost::algorithm::is_any_of(","),
                         boost::algorithm::token_compress_on);
            for (size_t j = 0; j < tokens.size(); j++) {
                (*deltas)[i][j] = boost::lexical_cast<int,string>(tokens[j]);
            }
            (*deltas)[i][tokens.size()] = -1;
        }

        // Print deltas
        cout << "Deltas to Encode: ";
        for (size_t i = 0; i < deltas_split.size(); i++) {
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
    }*/

/*static void ParseOptionXform(string str, int **xform_buf)
{
    vector<string> xform_split;

    boost::split(xform_split, str, boost::algorithm::is_any_of(","),
                 boost::algorithm::token_compress_on);

    *xform_buf = (int *) xmalloc(XFORM_MAX * sizeof(int));
    size_t index = 0;
    for (size_t i = 0; i < xform_split.size(); i++) {
        int option = XformOpt[xform_split[i]];
        if (option < 5) {
            (*xform_buf)[index++] = option;
            continue;
        }
        else if (option == 5) {  // BlockRow 
            for (size_t j = 7; j < 14; j++)
                (*xform_buf)[index++] = j;
        }
        else if (option == 6){  // BlockCol 
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
}
*/
