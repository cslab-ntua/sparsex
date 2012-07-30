/*
 * drle.cc -- Delta Run-Length Encoding Manager implementation
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <ios>
#include <exception>
#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>
#define FOREACH BOOST_FOREACH

#include "spm.h"
#include "drle.h"
#include "delta.h"

namespace bll = boost::lambda;

using namespace csx;

bool debug = false;

template<typename T> 
struct RLE {
    long freq;
    T val;
};

//
//  FIXME: MaxDelta() could be possibly merged with DeltaEncode()
// 
uint64_t MaxDelta(std::vector<uint64_t> xs)
{
    uint64_t prev_x, curr_x, max_delta;

    max_delta = 0;
    prev_x = curr_x = xs[0];
    for (size_t i = 1; i < xs.size(); ++i) {
        curr_x = xs[i];
        uint64_t delta = curr_x - prev_x;
        if (delta > max_delta)
            max_delta = delta;
        prev_x = curr_x;
    }

    return max_delta;
}

template<typename T>
T DeltaEncode(T input, size_t &max_delta)
{
    T output;
    typename T::iterator in, out;
    typename T::iterator::value_type prev, curr;

    output.resize(input.size());
    in = input.begin();
    out = output.begin();
    prev = *out++ = *in++;
    max_delta = 0;
    while (in < input.end()) {
        curr = *in++;
        size_t delta = curr - prev;
        *out++ = delta;
        if (delta > max_delta)
            max_delta = delta;
        prev = curr;
    }

    return output;
}

template<typename T>
T DeltaEncode(T input)
{
    size_t max_delta;
    return DeltaEncode(input, max_delta);
}

template<typename T>
std::vector<RLE<typename T::iterator::value_type> >
RLEncode(T input)
{
    typename T::iterator in;
    typename T::iterator::value_type curr;
    std::vector<RLE<typename T::iterator::value_type> > output;
    RLE<typename T::iterator::value_type> rle;

    in = input.begin();
    rle.freq = 1;
    rle.val = *in++;
    while (in < input.end()) {
        curr = *in;
        if (rle.val != curr) {
            output.push_back(rle);
            rle.freq = 1;
            rle.val = curr;
        } else {
            ++rle.freq;
        }

        ++in;
    }

    output.push_back(rle);
    return output;
}	                  

DRLE_Manager::DRLE_Manager(SPM *spm, long min_limit, long max_limit,
                           double min_perc, uint64_t sort_window_size,
                           split_alg_t split_type,
                           double sampling_portion,
                           uint64_t samples_max,
                           bool split_blocks,
                           bool onedim_blocks)
    : spm_(spm), min_limit_(min_limit),
      max_limit_(max_limit), min_perc_(min_perc),
      sort_window_size_(sort_window_size),
      split_type_(split_type),
      sampling_portion_(sampling_portion),
      samples_max_(samples_max),
      split_blocks_(split_blocks),
      onedim_blocks_(onedim_blocks)
{
    // These are delimiters, ignore them by default.
    AddIgnore(BLOCK_TYPE_START);
    AddIgnore(BLOCK_COL_START);
    AddIgnore(BLOCK_TYPE_END);
    
    if (!onedim_blocks_) {
        AddIgnore(BLOCK_ROW_TYPE_NAME(1));
        AddIgnore(BLOCK_COL_TYPE_NAME(1));
    }
    
    assert(sampling_portion_ >= 0 && sampling_portion_ <= 1 &&
           "invalid sampling portion");
    CheckAndSetSorting();
    if (sort_windows_) {
        if (sampling_portion > 0) {
            sort_window_size_ =
                sampling_portion_*spm_->GetNrNonzeros() / samples_max_;
        }
        
        ComputeSortSplits();
#if 0
        std::cout << "Window size (computed): "
                  << sort_window_size_ << std::endl;
        for (size_t i = 0; i < sort_splits_.size() - 1; ++i) {
            printf("(%ld, %ld)\n", sort_splits_[i+1] - sort_splits_[i],
                   sort_splits_nzeros_[i]);
        }
#endif

        if (samples_max_ > sort_splits_.size()) {
            samples_max_ = sort_splits_.size();
        }

        SelectSplits();
    } else {
        selected_splits_ = 0;
    }
    
    for (int i = 0; i < PRE_TIMER_END; i++)
        timer_init(&timers_[i]);
}
    
DeltaRLE::Stats DRLE_Manager::GenerateStats(uint64_t rs, uint64_t re)
{
    return GenerateStats(spm_, rs, re);
}

DeltaRLE::Stats DRLE_Manager::GenerateStats(SPM *spm, uint64_t rs, uint64_t re)
{
    std::vector<uint64_t> xs;
    DeltaRLE::Stats stats;
    std::vector<SpmRowElem*> elems;

    for (uint64_t i = rs; i < re; ++i) {
        for (SpmRowElem *elem = spm->RowBegin(i);
             elem != spm->RowEnd(i); ++elem) {
            if (elem->pattern == NULL) {
                xs.push_back(elem->x);
                elems.push_back(elem);
                continue;
            }

            UpdateStats(spm, xs, stats, elems);
        }

        UpdateStats(spm, xs, stats, elems);
    }

    return stats;
}

void DRLE_Manager::GenerateDeltaStats(SPM *spm, uint64_t rs, uint64_t re,
                                      DeltaRLE::Stats &stats)
{
    std::vector<uint64_t> xs;
    
    for (uint64_t i = rs; i < re; ++i) {
        int delta_bytes = 0;
        for (SpmRowElem *elem = spm->RowBegin(i);
             elem != spm->RowEnd(i); ++elem) {
            if (!elem->in_pattern)
                xs.push_back(elem->x);
            else
                elem->in_pattern = false;   // reset the marking

            if (elem->pattern_start) {
                // new pattern in the row creates a new delta unit
                if (xs.size()) {
                    delta_bytes = DeltaSize_getBytes(getDeltaSize(MaxDelta(xs)));
                    stats[0].npatterns +=
                        xs.size() / max_limit_ + (xs.size() % max_limit_ != 0);

                    stats[0].nnz += ((4 - delta_bytes)*xs.size()) / 4;
                    xs.clear();
                }
                elem->pattern_start = false;    // reset the marking
            }

            assert(!elem->pattern_start);
            assert(!elem->in_pattern);
        }
        
        if (xs.size()) {
            delta_bytes = DeltaSize_getBytes(getDeltaSize(MaxDelta(xs)));
            stats[0].npatterns +=
                xs.size() / max_limit_ + (xs.size() % max_limit_ != 0);
            stats[0].nnz += ((4 - delta_bytes)*xs.size()) / 4;
            xs.clear();
        }
    }
}

void DRLE_Manager::GenAllStats()
{
    DeltaRLE::Stats::iterator iter, tmp;
    DeltaRLE::Stats *sp;
    
    stats_.clear();
    deltas_to_encode_.clear();

    // first generate stats for the delta type
#ifdef SPM_HEUR_NEW
    GenerateDeltaStats(spm_, 0, spm_->GetNrRows(), stats_[NONE]);
#endif
    for (int t = HORIZONTAL; t != XFORM_MAX; ++t) {
        if (xforms_ignore_[t])
            continue;

        SpmIterOrder type = SpmTypes[t];

        // std::cout << "Checking for " << SpmTypesNames[t] << std::endl;
        if (sort_windows_ && spm_->GetNrRows() > samples_max_) {
            uint64_t samples_nnz = 0;
            for (size_t i = 0; i < samples_max_; ++i) {
                DeltaRLE::Stats w_stats;
                uint64_t window_start = sort_splits_[selected_splits_[i]];
                uint64_t window_size =
                    sort_splits_[selected_splits_[i]+1] - window_start;
                SPM *window = spm_->GetWindow(window_start, window_size);
                
                // std::cout << "wstart: " << window_start << std::endl;
                // std::cout << "wsize (rows): " << window_size << std::endl;
                // std::cout << "wsize (nnz): " << 
                //               sort_splits_nzeros_[selected_splits_[i]] <<
                //               std::endl;
                // Check for empty windows, since nonzeros might be captured
                // from previous patterns.
                if  (!window->nr_nzeros_)
                    goto exit_loop;
                    
                samples_nnz += sort_splits_nzeros_[selected_splits_[i]];
                window->Transform(type);
                w_stats = GenerateStats(window, 0, window->GetNrRows());
                UpdateStats(type, w_stats);
                window->Transform(HORIZONTAL);
                spm_->PutWindow(window);
            exit_loop:
                delete window;
            }
            
            CorrectStats(type, spm_->nr_nzeros_ / (double) samples_nnz);
        } else {
            spm_->Transform(type);
            stats_[type] = GenerateStats(0, spm_->GetNrRows());
            spm_->Transform(HORIZONTAL);
        }

        sp = &stats_[type];

        uint64_t block_align = IsBlockType(type);
        
        if (block_align && split_blocks_)
            CorrectBlockStats(sp, block_align);
        
        for (iter = sp->begin(); iter != sp->end(); ) {
            tmp = iter++;
            double p = (double) tmp->second.nnz / (double) spm_->nr_nzeros_;
            if (block_align == 1 && (unsigned int) tmp->first < min_limit_)
                // the dimension of one-dimensional blocks must exceed
                // min_limit_
                sp->erase(tmp);
            else if (p < min_perc_ || tmp->first >= CSX_PID_OFFSET)
                sp->erase(tmp);
            else
                deltas_to_encode_[type].insert(tmp->first);
        }

#ifdef SPM_HEUR_NEW
        GenerateDeltaStats(spm_, 0, spm_->GetNrRows(), *sp);
#endif
    }
}

namespace csx {

void DRLE_OutStats(DeltaRLE::Stats &stats, SPM &spm, std::ostream &os)
{
    DeltaRLE::Stats::iterator iter;

    for (iter = stats.begin(); iter != stats.end(); ++iter) {
        os << "    " << iter->first << "-> " << "np:"
           << iter->second.npatterns
           << " nnz: " << iter->second.nnz
           << " ("
           <<  100 * ((double) iter->second.nnz / (double) spm.GetNrNonzeros())
           << "%" << ")";
    }
}

} // end of csx namespace

void DRLE_Manager::AddIgnore(SpmIterOrder type)
{
    xforms_ignore_.set(type);
}

void DRLE_Manager::IgnoreAll()
{
    xforms_ignore_.set();
}

void DRLE_Manager::RemoveIgnore(SpmIterOrder type)
{
    // the following types are always ignored
    if (type <= NONE || type == BLOCK_TYPE_START ||
        (type == BLOCK_ROW_TYPE_NAME(1) && !onedim_blocks_) ||
        type == BLOCK_COL_START ||
        (type == BLOCK_COL_TYPE_NAME(1) && !onedim_blocks_) ||
        type == BLOCK_TYPE_END ||
        type >= XFORM_MAX)
        return;
    xforms_ignore_.reset(type);
}

void DRLE_Manager::RemoveAll()
{
    for (int t = NONE; t < XFORM_MAX; ++t)
        xforms_ignore_.reset(t);
}

SpmIterOrder DRLE_Manager::ChooseType()
{
    SpmIterOrder ret;
    long max_score;
    DRLE_Manager::StatsMap::iterator iter;

    ret = NONE;
    max_score = 0;
    for (iter = stats_.begin(); iter != stats_.end(); ++iter){
        long score = GetTypeScore(iter->first);
        if (score < 0) {
            AddIgnore(iter->first);
        } else if (score > max_score) {
            max_score = score;
            ret = iter->first;
        }
    }

    return ret;
}

long DRLE_Manager::GetTypeScore(SpmIterOrder type)
{
    DeltaRLE::Stats *sp;
    DeltaRLE::Stats::iterator iter;
    long ret;

    ret = 0;
    if (stats_.find(type) == stats_.end())
        return ret;
        
    sp = &stats_[type];

    long nr_nzeros_encoded = 0;
    long nr_patterns = 0;
    long nr_deltas = spm_->nr_deltas_;

    for (iter = sp->begin(); iter != sp->end(); ++iter) {
        if (iter->first == 0 && type != NONE) {
            nr_deltas += iter->second.npatterns;
        } else {
            nr_nzeros_encoded += iter->second.nnz;
            nr_patterns += iter->second.npatterns;
        }
    }

#ifdef SPM_HEUR_NEW
    long nr_switches;
    if (type == NONE)
        nr_switches = nr_deltas;
    else
        nr_switches = nr_deltas + nr_patterns;

    double coeff;
    char *coeff_str = getenv("_COEFF");
    if (coeff_str)
        coeff = atof(coeff_str);
    else
        coeff = 1.;

     ret = (nr_nzeros_encoded - nr_patterns)*coeff - nr_switches;
#else
     ret = nr_nzeros_encoded - nr_patterns;
#endif
     return ret;
}

void DRLE_Manager::Encode(SpmIterOrder type)
{
    SPM::Builder *SpmBld;
    SpmIterOrder oldtype;
    std::vector<SpmRowElem> new_row;
    uint64_t nr_size;
    SpmRowElem *elems;

    if (type == NONE && ((type = ChooseType()) == NONE))
        return;

    SpmBld = new SPM::Builder(spm_, spm_->elems_size_, spm_->max_rowptr_size_);

    // Transform matrix to the desired iteration order
    oldtype = spm_->type_;
    spm_->Transform(type);

    for (uint64_t i = 0; i < spm_->GetNrRows(); ++i) {
        EncodeRow(spm_->RowBegin(i), spm_->RowEnd(i), new_row);
        nr_size = new_row.size();
        if (nr_size > 0) {
            timer_start(&timers_[PRE_TIMER_ALLOC]);
            elems = SpmBld->AllocElems(nr_size);
            timer_pause(&timers_[PRE_TIMER_ALLOC]);
            for (uint64_t i = 0; i < nr_size; ++i)
                MakeRowElem(new_row[i], elems + i);
        }

        new_row.clear();
        SpmBld->NewRow();
    }

    SpmBld->Finalize();
    delete SpmBld;
    
    // Transform matrix to the original iteration order
    spm_->Transform(oldtype);
    AddIgnore(type);
}

void DRLE_Manager::Decode(SpmIterOrder type)
{
    SPM::Builder *SpmBld;
    SpmIterOrder oldtype;
    std::vector<SpmRowElem> new_row;
    uint64_t nr_size;
    SpmRowElem *elems;

    if (type == NONE)
        return;

    // Transform matrix to the desired iteration order
    oldtype = spm_->type_;
    spm_->Transform(type);

    // Do the decoding
    SpmBld = new SPM::Builder(spm_, spm_->elems_size_, spm_->max_rowptr_size_);
    for (uint64_t i = 0; i < spm_->GetNrRows(); ++i) {
        DecodeRow(spm_->RowBegin(i), spm_->RowEnd(i), new_row);
        nr_size = new_row.size();
        if (nr_size > 0) {
            elems = SpmBld->AllocElems(nr_size);
            for (uint64_t i = 0; i < nr_size; ++i)
                MakeRowElem(new_row[i], elems + i);
	}

	new_row.clear();
        SpmBld->NewRow();
    }

    SpmBld->Finalize();
    delete SpmBld;

    // Transform matrix to the original iteration order
    spm_->Transform(oldtype);
    RemoveIgnore(type);
}

void DRLE_Manager::EncodeAll(std::ostream &os)
{
    SpmIterOrder type = NONE;
    StatsMap::iterator iter; 
    SpmIterOrder enc_seq[22];
    int counter = 0;
    
    timer_start(&timers_[PRE_TIMER_TOTAL]);
    for (;;) {
        timer_start(&timers_[PRE_TIMER_STATS]);
        GenAllStats();
        timer_pause(&timers_[PRE_TIMER_STATS]);
        OutStats(os);
        type = ChooseType();
        if (type == NONE)
            break;
        timer_start(&timers_[PRE_TIMER_ENCODE]);
        os << "Encode to " << SpmTypesNames[type] << std::endl;
        
        Encode(type);

        // Add the deltas produced by this encoding
        spm_->nr_deltas_ += stats_[type][0].npatterns;

        enc_seq[counter++] = type;
        timer_pause(&timers_[PRE_TIMER_ENCODE]);
    }

    timer_pause(&timers_[PRE_TIMER_TOTAL]);

    os << "Encoding sequence: ";
    if (counter == 0)
        os << "NONE";
    else
        os << SpmTypesNames[enc_seq[0]];

    for (int i = 1; i < counter; i++)
        os << ", " << SpmTypesNames[enc_seq[i]];
    os << std::endl;
    
    /*
    double t;

    for (int i = 0; i < PRE_TIMER_END; i++) {
        t = timer_secs(&timers_[i]);
        os << timers_desc_[i] << t << std::endl;
    }
    */
}

void DRLE_Manager::MakeEncodeTree()
{
    uint32_t i,j;
    uint32_t count = 0;
    StatsMap::iterator iter;
    std::vector<Node> nodes;
    std::vector<Node>::iterator it,it2;

    nodes.push_back(Node(0));
    while (nodes.size() != 0) {
        for (i = 0; nodes[0].type_ignore_[i] != NONE; ++i)
            AddIgnore(nodes[0].type_ignore_[i]);
            
        GenAllStats();
        for (i = 0; nodes[0].type_ignore_[i] != NONE; ++i)
            RemoveIgnore(nodes[0].type_ignore_[i]);
            
        for (iter = stats_.begin(); iter != stats_.end(); ++iter)
            if (GetTypeScore(iter->first) == 0)
                nodes[0].Ignore(iter->first);
                
        it = nodes.begin();
        for (iter = stats_.begin(); iter != stats_.end(); ++iter) {
            if (GetTypeScore(iter->first) > 0) {
                ++it;
                it = nodes.insert(it, nodes[0].MakeChild(iter->first,
                                  deltas_to_encode_[iter->first]));
            }
        }

        if (it != nodes.begin()) {
            Encode(nodes[1].type_path_[nodes[1].depth_-1]);
        } else {
            nodes[0].PrintNode();
            ++count;
            if (nodes.size() != 1) {
                i = 0;
                while (nodes[0].type_path_[i] == nodes[1].type_path_[i])
                    ++i;
                for (j = 1; j <= nodes[0].depth_ - i; ++j)
                     Decode(nodes[0].type_path_[nodes[0].depth_-j]);
                     
                while (i < nodes[1].depth_) {
                    Encode(nodes[1].type_path_[i]);
                    ++i;
                }
            } else {
                for (j = 1; j <= nodes[0].depth_; ++j)
                    Decode(nodes[0].type_path_[nodes[0].depth_-j]);
            }
        }

        delete nodes[0].type_path_;
        delete nodes[0].type_ignore_;
        nodes.erase(nodes.begin());
    }

    std::cout << "Tree has " << count << " possible paths" << std::endl;
    exit(0);
}

void DRLE_Manager::EncodeSerial(int *xform_buf, int **deltas)
{
    for (uint32_t i = 0; i < XFORM_MAX; ++i)
        AddIgnore((SpmIterOrder) i);

    for (int i = 0; xform_buf[i] != -1; ++i) {
        SpmIterOrder t = static_cast<SpmIterOrder>(xform_buf[i]);

        RemoveIgnore(t);
        for (int j = 0; deltas[i][j] != -1; ++j)
            deltas_to_encode_[t].insert(deltas[i][j]);
            
        Encode(t);
        AddIgnore(t);
    }
}

void DRLE_Manager::OutputSortSplits(std::ostream& out)
{
    sort_split_iterator iter;
    for (iter = sort_splits_.begin(); iter != sort_splits_.end() - 1; ++iter) {
        uint64_t rs = *iter;
        uint64_t re = *(iter + 1);
        uint64_t nnz = spm_->rowptr_[re] - spm_->rowptr_[rs];

        out << "(rs, re, nnz) = (" << rs << ", " << re << ", " << nnz << ")"
            << std::endl;
    }
}

void DRLE_Manager::DoEncode(std::vector<uint64_t> &xs, std::vector<double> &vs,
                            std::vector<SpmRowElem> &encoded)
{
    uint64_t col;
    std::vector< RLE<uint64_t> > rles;
    std::set<uint64_t> *deltas_set;
    std::vector<double>::iterator vi = vs.begin();
    SpmRowElem elem;

    if (IsBlockType(spm_->type_)) {
        if (!split_blocks_)
            DoEncodeBlock(xs, vs, encoded);
        else
            DoEncodeBlockAlt(xs, vs, encoded);
        return;
    }

    // do a delta run-length encoding of the x values
    rles = RLEncode(DeltaEncode(xs));

    // Not all delta rles are to be encoded, only those
    // that are in the ->deltas_to_encode set
    deltas_set = &deltas_to_encode_[spm_->type_];

    col = 0;
    elem.pattern = NULL;
    FOREACH (RLE<uint64_t> rle, rles) {
        long rle_freq, rle_start;
        rle_freq = rle.freq;
        if (rle_freq != 1 && deltas_set->find(rle.val) != deltas_set->end()) {

            col += rle.val;
            if (col != rle.val) {
                SpmRowElem *last_elem = &encoded.back();
                rle_start = col;
                rle_freq = rle.freq;
#ifndef SPM_NUMA
                if (last_elem->pattern == NULL) {
                    // include the previous element, too
                    rle_start -= rle.val;
                    rle_freq++;
                    encoded.pop_back();
                    --vi;
                }
#endif
            } else {
                // we are the first unit in the row
                rle_start = col;
                rle_freq = rle.freq;
            }
            
            int nr_parts = 0;
            while (rle_freq >= min_limit_) {
                SpmRowElem *last_elem;
                long curr_freq = std::min(max_limit_, rle_freq);
                
                elem.x = rle_start;
                encoded.push_back(elem);
                last_elem = &encoded.back();
                last_elem->pattern = new DeltaRLE(curr_freq, rle.val,
                                                  spm_->type_);
                last_elem->vals = new double[curr_freq];
                std::copy(vi, vi + curr_freq, last_elem->vals);
                vi += curr_freq;
                rle_start += rle.val * curr_freq;
                rle_freq -= curr_freq;
                ++nr_parts;
            }

            // leave col at the last element of the pattern
            col = rle_start - rle.val;
        }

        // add individual elements
        for (int i = 0; i < rle_freq; ++i) {
            col += rle.val;
            elem.x = col;
            elem.val = *vi++;
            encoded.push_back(elem);
        }
    }

    assert(vi == vs.end() && "not all elements processed or out of bounds");
    xs.clear();
    vs.clear();
}

//static uint64_t nr_lines = 0;
void DRLE_Manager::DoEncodeBlock(std::vector<uint64_t> &xs,
                                 std::vector<double> &vs,
                                 std::vector<SpmRowElem> &encoded)
{
    uint64_t col;
    std::vector< RLE<uint64_t> > rles;
    std::set<uint64_t> *deltas_set;
    std::vector<double>::iterator vi = vs.begin();
    SpmRowElem elem;
    SpmRowElem *last_elem;

    // do a delta run-length encoding of the x values
    rles = RLEncode(DeltaEncode(xs));

    // Not all delta rles are to be encoded, only those
    // that are in the ->deltas_to_encode set
    deltas_set = &deltas_to_encode_[spm_->type_];

    int block_align = IsBlockType(spm_->type_);
    assert(block_align && "not a block type");

    col = 0;
    elem.pattern = NULL;
    FOREACH (RLE<uint64_t> rle, rles) {
        uint64_t skip_front, skip_back, nr_elem;

        col += rle.val;
        if (col == 1) {
            skip_front = 0;
            nr_elem = rle.freq;
        } else {
            skip_front = (col - 2) % block_align;
            if (skip_front != 0)
                skip_front  = block_align - skip_front;
                
            nr_elem = rle.freq + 1;
        }

        if (nr_elem > skip_front)
            nr_elem -= skip_front;
        else
            nr_elem = 0;
            
        skip_back = nr_elem % block_align;
        if (nr_elem > skip_back)
            nr_elem -= skip_back;
        else
            nr_elem = 0;
            
        if (rle.val == 1 &&
            deltas_set->find(nr_elem / block_align) != deltas_set->end() &&
            nr_elem >= (uint64_t) 2 * block_align) {

            uint64_t rle_start;

            if (col != 1) {
                rle_start = col - 1;
                encoded.pop_back();
                --vi;
            } else {
                rle_start = col;
            }

            // Add elements skipped from start
            for (uint64_t i = 0; i < skip_front; ++i) {
                elem.x = rle_start + i;
                elem.val = *vi++;
                encoded.push_back(elem);
            }

            // Align max_limit
            uint64_t max_limit = max_limit_ / block_align * block_align;
            uint64_t nr_blocks = nr_elem / max_limit;
            uint64_t nr_elem_block = std::min(max_limit, nr_elem);

            if (nr_blocks == 0)
                nr_blocks = 1;
            else
                // Adjust skip_back with leftover items from max_limit alignment
                skip_back += nr_elem - nr_elem_block * nr_blocks;

            for (uint64_t i = 0; i < nr_blocks; ++i) {
                // Add the blocks
                elem.x = rle_start + skip_front + i * nr_elem_block;
                encoded.push_back(elem);
                last_elem = &encoded.back();
                last_elem->pattern = new BlockRLE(nr_elem_block,
                                                  nr_elem_block / block_align,
                                                  spm_->type_);
                last_elem->vals = new double[nr_elem_block];
                std::copy(vi, vi + nr_elem_block, last_elem->vals);
                vi += nr_elem_block;
            }

            // Add the remaining elements
            for (uint64_t i = 0; i < skip_back; ++i) {
                elem.x = rle_start + skip_front + nr_elem_block * nr_blocks + i;
                elem.val = *vi++;
                encoded.push_back(elem);
            }
        } else {
            // add individual elements
            for (int i = 0; i < rle.freq; ++i) {
                elem.x = col + i * rle.val;
                elem.val = *vi++;
                encoded.push_back(elem);
            }
        }

        col += rle.val * (rle.freq - 1);
    }

    assert(vi == vs.end() && "out of bounds");
    xs.clear();
    vs.clear();
}

void DRLE_Manager::DoEncodeBlockAlt(std::vector<uint64_t> &xs,
                                    std::vector<double> &vs,
                                    std::vector<SpmRowElem> &encoded)
{
    uint64_t col;
    std::vector< RLE<uint64_t> > rles;
    std::set<uint64_t> *deltas_set;
    std::vector<double>::iterator vi = vs.begin();
    SpmRowElem elem;
    SpmRowElem *last_elem;

    // do a delta run-length encoding of the x values
    rles = RLEncode(DeltaEncode(xs));

    // Not all delta rles are to be encoded, only those
    // that are in the ->deltas_to_encode set
    deltas_set = &deltas_to_encode_[spm_->type_];
    int block_align = IsBlockType(spm_->type_);
    assert(block_align && "not a block type");

    col = 0;
    elem.pattern = NULL;
    FOREACH (RLE<uint64_t> rle, rles) {
        uint64_t skip_front, skip_back, nr_elem;

        col += rle.val;
        if (col == 1) {
            skip_front = 0;
            nr_elem = rle.freq;
        } else {
            skip_front = (col - 2) % block_align;
            if (skip_front != 0)
                skip_front  = block_align - skip_front;
            nr_elem = rle.freq + 1;
        }

        if (nr_elem > skip_front)
            nr_elem -= skip_front;
        else
            nr_elem = 0;
        skip_back = nr_elem % block_align;
        nr_elem -= skip_back;
        if (rle.val == 1 && nr_elem >= (uint64_t) 2 * block_align) {
            uint64_t rle_start;

            if (col != 1) {
                rle_start = col - 1;
                encoded.pop_back();
                --vi;
            } else {
                rle_start = col;
            }

            // Add elements skipped from start
            for (uint64_t i = 0; i < skip_front; ++i) {
                elem.x = rle_start++;
                elem.val = *vi++;
                encoded.push_back(elem);
            }

            // Split Blocks to Match the Appropriate Sizes for Encode
            uint64_t other_dim = nr_elem / block_align;

            for (std::set<uint64_t>::reverse_iterator i = deltas_set->rbegin();
                                              i != deltas_set->rend(); ++i) {
                while (other_dim >= *i) {
                    //We have a new block RLE
                    uint64_t nr_elem_block = block_align * (*i);
                    
                    elem.x = rle_start;
                    rle_start += nr_elem_block;
                    encoded.push_back(elem);
                    last_elem = &encoded.back();
                    last_elem->pattern = new BlockRLE(nr_elem_block, *i,
                                                      spm_->type_);
                    last_elem->vals = new double[nr_elem_block];
                    std::copy(vi, vi + nr_elem_block, last_elem->vals);
                    vi += nr_elem_block;
                    nr_elem -= nr_elem_block;
                    other_dim -= (*i);
                }
            }

            // Add the remaining elements
            skip_back += nr_elem;
            for (uint64_t i = 0; i < skip_back; ++i) {
                elem.x = rle_start++;
                elem.val = *vi++;
                encoded.push_back(elem);
            }
        } else {
            // add individual elements
            for (int i = 0; i < rle.freq; ++i) {
                elem.x = col + i * rle.val;
                elem.val = *vi++;
                encoded.push_back(elem);
            }
        }

        col += rle.val * (rle.freq - 1);
    }

    assert(vi == vs.end() && "out of bounds");
    xs.clear();
    vs.clear();
}

void DRLE_Manager::DoDecode(const SpmRowElem *elem, 
                            std::vector<SpmRowElem> &newrow)
{
    long int i;
    uint64_t cur_x;
    SpmRowElem new_elem;

    new_elem.pattern = NULL;
    cur_x = elem->x;
    for (i = 0; i < elem->pattern->GetSize(); ++i) {
        new_elem.x = cur_x;
        new_elem.val = elem->vals[i];
        newrow.push_back(new_elem);
        cur_x = elem->pattern->GetNextCol(cur_x);
    }

    delete elem->pattern;
    delete elem->vals;
}

void DRLE_Manager::EncodeRow(const SpmRowElem *rstart, const SpmRowElem *rend,
                             std::vector<SpmRowElem> &newrow)
{
    std::vector<uint64_t> xs;
    std::vector<double> vs;

    // gather x values into xs vector until a pattern is found
    // and encode them using DoEncode()
    for (const SpmRowElem *e = rstart; e < rend; ++e) {
        if (e->pattern == NULL) {
            xs.push_back(e->x);
            vs.push_back(e->val);
            continue;
        }

        if (xs.size() != 0)
            DoEncode(xs, vs, newrow);

        newrow.push_back(*e);
        
        if (e->pattern)
            delete e->pattern;
    }

    // Encode any remaining elements
    if (xs.size() != 0)
        DoEncode(xs, vs, newrow);
}

void DRLE_Manager::DecodeRow(const SpmRowElem *rstart,
                             const SpmRowElem *rend,
                             std::vector<SpmRowElem> &newrow)
{
    for (const SpmRowElem *e = rstart; e < rend; ++e) {
        if (e->pattern != NULL && e->pattern->GetType() == spm_->type_) {
            DoDecode(e, newrow);
        } else {
            newrow.push_back(*e);
            if (e->pattern)
                delete e->pattern;
        }
    }
}

void DRLE_Manager::CorrectBlockStats(DeltaRLE::Stats *sp, uint64_t block_align)
{
    DeltaRLE::Stats temp, temp2;
    uint64_t max_block_dim;
    DeltaRLE::Stats::iterator tmp;
    DeltaRLE::Stats::reverse_iterator iter, rtmp;

    max_block_dim = max_limit_ / block_align;
    temp[max_block_dim].nnz = 0;
    temp[max_block_dim].npatterns = 0;
    for (iter = sp->rbegin(); iter != sp->rend() &&
         iter->first * block_align > (unsigned) max_limit_; ++iter) {
        
        uint64_t div_factor = iter->first / max_block_dim;
        uint64_t mod_factor = iter->first % max_block_dim;
        uint64_t max_block_size = max_block_dim * block_align;
        uint64_t block_patterns = div_factor * iter->second.npatterns;
        uint64_t remaining_nnz;
        temp[max_block_dim].nnz += block_patterns * max_block_size;
        temp[max_block_dim].npatterns += block_patterns;
        remaining_nnz = iter->second.nnz - block_patterns * max_block_size;
        if (mod_factor >= 2) {
            tmp = temp.find(mod_factor);
            if (tmp == temp.end()) {
                temp[mod_factor].nnz = remaining_nnz;
                temp[mod_factor].npatterns = iter->second.npatterns;
            } else {
                temp[mod_factor].nnz += remaining_nnz;
                temp[mod_factor].npatterns += iter->second.npatterns;
            }
        }
    }
    
    if (temp[max_block_dim].npatterns == 0) {
        assert(temp[max_block_dim].nnz == 0);
        temp.erase(max_block_dim);
    }
    
    for ( ; iter != sp->rend(); ++iter) {
        tmp = temp.find(iter->first);
        if (tmp == temp.end()) {
            temp[iter->first].nnz = iter->second.nnz;
            temp[iter->first].npatterns = iter->second.npatterns;
        } else {
            temp[iter->first].nnz += iter->second.nnz;
            temp[iter->first].npatterns += iter->second.npatterns;
        }
    }

    uint64_t div_factor, mod_factor;
    
    rtmp = sp->rbegin();
    for (iter = temp.rbegin(); iter != temp.rend(); ++iter) {
        if (((double) iter->second.nnz) / ((double) spm_->nr_nzeros_) >= 
            min_perc_) {
            for ( ; rtmp != sp->rend() && rtmp->first >= iter->first; ++rtmp) {
                div_factor = rtmp->first / iter->first;
                mod_factor = rtmp->first % iter->first;
               
                tmp = temp2.find(iter->first);
                if (tmp == temp2.end()) {
                    temp2[iter->first].npatterns =
                        div_factor * rtmp->second.npatterns;
                    temp2[iter->first].nnz =
                        div_factor * (rtmp->second.nnz / rtmp->first) *
                        iter->first;
                } else {   
                    tmp->second.npatterns +=
                        div_factor * rtmp->second.npatterns;
                    tmp->second.nnz +=
                        div_factor * (rtmp->second.nnz / rtmp->first) *
                        iter->first;
                }

                if (mod_factor >= 2) {
                    tmp = temp2.find(mod_factor);
                    if (tmp == temp2.end()) {
                        temp2[mod_factor].npatterns = rtmp->second.npatterns;
                        temp2[mod_factor].nnz =
                            (rtmp->second.nnz / rtmp->first) * mod_factor;
                    } else {
                        tmp->second.npatterns += rtmp->second.npatterns;
                        tmp->second.nnz +=
                            (rtmp->second.nnz / rtmp->first) * mod_factor;
                    }
                }
            }
        }
    }
    sp->clear();
    *sp = temp2;
}

void DRLE_Manager::HandleStats(DeltaRLE::Stats *sp, uint64_t size,
                               uint64_t block_align)
{
    DeltaRLE::Stats temp;
    DeltaRLE::Stats::iterator iter, not_passed, tmp, tmp2;

    not_passed = sp->end();
    for (iter = sp->end(); iter != sp->begin(); ) {
        double p = 0;

        iter--;
        if (iter->first * block_align <= (unsigned) max_limit_) {
            for (tmp = iter; tmp != not_passed; ++tmp)
                p += ((double) tmp->second.npatterns * block_align * 
                     iter->first * (tmp->first / iter->first) / (double) size);
                     
            if (p > min_perc_) {
                uint64_t block_size = iter->first * block_align;

                temp[iter->first].nnz = iter->second.nnz;
                temp[iter->first].npatterns = iter->second.npatterns;
                tmp = iter;
                ++tmp;
                for (; tmp != not_passed; ++tmp) {
                    uint64_t div_factor = tmp->first / iter->first;
                    uint64_t mod_factor = tmp->first % iter->first;
                    temp[iter->first].nnz +=
                        tmp->second.npatterns * div_factor * block_size;
                    temp[iter->first].npatterns +=
                        tmp->second.npatterns * div_factor;
                    tmp->second.nnz -=
                        tmp->second.npatterns * div_factor * block_size;
                    if (mod_factor >= 2) {
                        tmp2 = temp.find(mod_factor);
                        if (tmp2 == temp.end()) {
                            temp[mod_factor].nnz = tmp->second.nnz;
                            temp[mod_factor].npatterns = tmp->second.npatterns;
                        } else {
                            temp[mod_factor].nnz = tmp->second.nnz;
                            temp[mod_factor].npatterns = tmp->second.npatterns;
                        }
                    }
                }

                not_passed = iter;
            }
        }
    }

    for (iter = sp->begin(); iter != sp->end(); ++iter)
        sp->erase(iter);
        
    for (iter = temp.begin(); iter != temp.end(); ++iter)
        sp->insert(std::pair<const long unsigned int,
                   csx::DeltaRLE::StatsVal>(iter->first,iter->second));
}

void DRLE_Manager::UpdateStats(SpmIterOrder type, DeltaRLE::Stats stats)
{
    if (stats_.find(type) == stats_.end()) {
        stats_[type] = stats;
    } else {
        for (DeltaRLE::Stats::const_iterator it = stats.begin();
            it != stats.end(); ++it) {
            stats_[type][it->first].Update(it->second);
        }
    }
}

void DRLE_Manager::UpdateStats(SPM *spm, std::vector<uint64_t> &xs,
                               DeltaRLE::Stats &stats,
                               std::vector<SpmRowElem*> &elems)
{
    std::vector< RLE<uint64_t> > rles;
    uint64_t block_align = IsBlockType(spm->type_);

    if (block_align) {
        DRLE_Manager::UpdateStatsBlock(xs, stats, block_align, elems);
        return;
    }

    assert(xs.size() == elems.size());
    if (xs.size() == 0)
        return;

    rles = RLEncode(DeltaEncode(xs));
    uint64_t col = 0;
    bool last_rle_patt = false; // turn on for encoded units
    std::vector<SpmRowElem*>::iterator ei = elems.begin();
    std::vector<SpmRowElem*>::iterator ee = elems.end();
    FOREACH(RLE<uint64_t> &rle, rles) {
#ifndef SPM_NUMA
        long real_limit = (col && !last_rle_patt) ? min_limit_ - 1 : min_limit_;
#else
        long real_limit = min_limit_;
#endif
        if (rle.freq > 1 && rle.freq >= real_limit) {
#ifndef SPM_NUMA
            long real_nnz = (col && !last_rle_patt) ?
                rle.freq + 1 : rle.freq;
#else
            long real_nnz = rle.freq;
#endif
            long rem_nnz = real_nnz % max_limit_;
            stats[rle.val].nnz += real_nnz;
            stats[rle.val].npatterns += real_nnz / max_limit_ + (rem_nnz != 0);
            if (rem_nnz && rem_nnz < min_limit_) {
                // remove the pattern for the remainder elements; too short
                stats[rle.val].npatterns--;
                stats[rle.val].nnz -= rem_nnz;
            } else {
                rem_nnz = 0;
            }

            if (real_nnz > rle.freq)
                --ei;

            // Mark the elements of the pattern and its start
            (*ei)->pattern_start = true;
            for (long i = 0; i < real_nnz - rem_nnz; ++i) {
                assert(ei < ee);
                (*ei)->in_pattern = true;
                ++ei;
            }

            last_rle_patt = true;
        } else {
            last_rle_patt = false;
            ei += rle.freq;
        }

        col += rle.val;
    }

    xs.clear();
    elems.clear();
}

void DRLE_Manager::UpdateStatsBlock(std::vector<uint64_t> &xs,
                                    DeltaRLE::Stats &stats,
                                    uint64_t block_align,
                                    std::vector<SpmRowElem*> &elems)
{
    std::vector<RLE<uint64_t> > rles;

    assert(block_align && "not a block type");
    if (xs.size() == 0)
        return;
        
    rles = RLEncode(DeltaEncode(xs));

    uint64_t unit_start = 0;
    std::vector<SpmRowElem*>::iterator ei = elems.begin();
    std::vector<SpmRowElem*>::iterator ee = elems.end();
    FOREACH (RLE<uint64_t> &rle, rles) {
        unit_start += rle.val;
        if (rle.val == 1) {
            // Start of the real block is at `unit_start - 1' with
            // one-based indexing. When computing the `%' we need
            // zero-based indexing.
            long nr_elem;
            long skip_front;

            if (unit_start == 1) {
                skip_front = 0;
                nr_elem = rle.freq;
            } else {
                skip_front = (unit_start-2) % block_align;
                if (skip_front != 0)
                    skip_front = block_align - skip_front;
                nr_elem = rle.freq + 1;
            }
            
            if (nr_elem > rle.freq)
                --ei;

            if (nr_elem > skip_front)
                nr_elem -= skip_front;
            else
                nr_elem = 0;

            long other_dim = nr_elem / block_align;

            ei += std::min(skip_front, nr_elem);
            if (other_dim >= 2) {
                stats[other_dim].nnz += other_dim * block_align;
                stats[other_dim].npatterns++;

                // Mark the elements
                (*ei)->pattern_start = true;
                for (uint64_t i = 0; i < nr_elem; ++i) {
                    assert(ei < ee);
                    if (i < stats[other_dim].nnz)
                        (*ei)->in_pattern = true;
                    ++ei;
                }
            }
        } else {
            ei += rle.freq;
        }
        
        unit_start += rle.val * (rle.freq - 1);
    }

    xs.clear();
    elems.clear();
}

void DRLE_Manager::CorrectStats(SpmIterOrder type, double factor)
{
    if (stats_.find(type) == stats_.end())
        return;

    DeltaRLE::Stats *l_stats = &stats_[type];
    DeltaRLE::Stats::iterator iter;

    for (iter = l_stats->begin(); iter != l_stats->end(); ++iter) {
        iter->second.nnz = (uint64_t) (iter->second.nnz * factor);
        iter->second.npatterns = (long) (iter->second.npatterns * factor);
        assert(iter->second.nnz <= spm_->GetNrNonzeros() &&
               "nonzeros of pattern exceed nonzeros of matrix");
    }
}

void DRLE_Manager::SelectSplits()
{
    size_t nr_splits = sort_splits_.size();
    size_t nr_samples = samples_max_;

    selected_splits_ = new size_t[nr_samples];
    if (nr_samples == nr_splits) {
        for (size_t i = 0; i < nr_splits; ++i)
            selected_splits_[i] = i;
        return;
    }

    if (nr_samples > nr_splits / 2) {
        for (size_t i = 0; i < nr_splits/2; ++i)
            selected_splits_[i] = i;
        nr_samples -= nr_splits / 2;
        nr_splits -= nr_splits / 2;
    }

    size_t skip = nr_splits / (nr_samples+1);
    for (size_t i = 0; i < nr_samples; ++i) {
        selected_splits_[i] = (i+1)*skip;
    }
    
    return;
}


void DRLE_Manager::ComputeSortSplits()
{
    switch (split_type_) {
    case DRLE_Manager::SPLIT_BY_ROWS:
        DoComputeSortSplitsByRows();
        break;
    case DRLE_Manager::SPLIT_BY_NNZ:
        DoComputeSortSplitsByNNZ();
        break;
    default:
        assert(false && "unknown split algorithm");
    }
}

void DRLE_Manager::CheckAndSetSorting()
{
    switch (split_type_) {
    case DRLE_Manager::SPLIT_BY_ROWS:
        DoCheckSortByRows();
        break;
    case DRLE_Manager::SPLIT_BY_NNZ:
        DoCheckSortByNNZ();
        break;
    default:
        assert(false && "unknown split algorithm");
    }
}

void DRLE_Manager::DoComputeSortSplitsByRows()
{
    uint64_t nr_rows = spm_->GetNrRows();
    uint64_t i;

    for (i = 0; i <= nr_rows; i += sort_window_size_)
        sort_splits_.push_back(i);
        
    if (i > nr_rows && i - nr_rows < sort_window_size_ / 2) {
        sort_splits_.push_back(nr_rows);
    } else {
        sort_splits_.pop_back();
        sort_splits_.push_back(nr_rows);
    }
}

void DRLE_Manager::DoComputeSortSplitsByNNZ()
{
    uint64_t nzeros_cnt;
    uint64_t nr_rows = spm_->GetNrRows();

    nzeros_cnt = 0;
    sort_splits_.push_back(0);
    for (uint64_t i = 0; i < nr_rows; ++i) {
        uint64_t new_nzeros_cnt = 
            nzeros_cnt + spm_->rowptr_[i+1] - spm_->rowptr_[i];
        if (new_nzeros_cnt < sort_window_size_) {
            nzeros_cnt = new_nzeros_cnt;
        } else {
            sort_splits_.push_back(i + 1);
            sort_splits_nzeros_.push_back(new_nzeros_cnt);
            nzeros_cnt = 0;
        }
    }

    if (nzeros_cnt) {
        uint64_t &last_split_nz = sort_splits_nzeros_.back();
        last_split_nz += nzeros_cnt;
        if (nzeros_cnt > sort_window_size_ / 2) {
            sort_splits_.push_back(nr_rows);
        } else {
            sort_splits_.pop_back();
            sort_splits_.push_back(nr_rows);
            
        }    
    }
}

void DRLE_Manager::DoCheckSortByRows()
{
    assert(sort_window_size_ <= spm_->GetNrRows() && "invalid sort window");
    if (sampling_portion_ == 0.0 || sort_window_size_ == spm_->GetNrRows())
        sort_windows_ = false;
    else
        sort_windows_ = true;
}

void DRLE_Manager::DoCheckSortByNNZ()
{
    assert(sort_window_size_ <= spm_->elems_size_ && "invalid sort window");
    if (sampling_portion_ == 0 || sort_window_size_ == spm_->elems_size_)
        sort_windows_ = false;
    else
        sort_windows_ = true;
}


void DRLE_Manager::OutStats(std::ostream &os)
{
    DRLE_Manager::StatsMap::iterator iter;
    for (iter = stats_.begin(); iter != stats_.end(); ++iter){
#if SPM_HEUR_NEW
        bool valid_stat = (iter->first == 0) || (iter->first != 0 && iter->second.size() > 1);
#else
        bool valid_stat = !iter->second.empty();
#endif
        if (valid_stat) {
             os << SpmTypesNames[iter->first] << "\ts:" 
                << GetTypeScore(iter->first);
             DRLE_OutStats(iter->second, *(spm_), os);
             os << std::endl;
         }
    }
}

Node::Node(uint32_t depth) : depth_(depth)
{
    uint32_t i;
    
    type_path_ = new SpmIterOrder[depth_];
    type_ignore_ = new SpmIterOrder[XFORM_MAX];
    for (i = 0; i < ((uint32_t) XFORM_MAX); i++)
        type_ignore_[i] = NONE;
}

void Node::PrintNode()
{
    for (uint32_t i = 0; i < depth_; ++i) {
        if (i != 0)
            std::cout << ",";
            
        std::cout << type_path_[i];
    }

    std::cout << std::endl;
    for (uint32_t i = 0; i < depth_; ++i) {
        SpmIterOrder temp_type = type_path_[i];
        std::set<uint64_t>::iterator it = deltas_path_[temp_type].begin();

        if (i != 0)
            std::cout << ",";
            
        std::cout << "{";
        for (uint32_t i = 1;
             i < static_cast<uint32_t>(deltas_path_[temp_type].size());
             ++i) {
            std::cout << *it << ",";
            ++it;
        }

        std::cout << *it << "}";
    }

    std::cout << std::endl;
}

void Node::Ignore(SpmIterOrder type)
{
    uint32_t i = 0;
    
    while (type_ignore_[i] != NONE) {
        assert(type_ignore_[i] != type && "type already ignored");
        ++i;
    }
    
    type_ignore_[i] = type;
}

Node Node::MakeChild(SpmIterOrder type, std::set<uint64_t> deltas)
{
    Node new_node = Node(depth_ + 1);

    for (uint32_t i = 0; i < static_cast<uint32_t>(XFORM_MAX); ++i)
        new_node.type_ignore_[i] = type_ignore_[i];

    for (uint32_t i = 0; i < depth_; ++i) {
        SpmIterOrder temp_type = type_path_[i];
        new_node.type_path_[i] = temp_type;
        new_node.deltas_path_[temp_type] = deltas_path_[temp_type];
    }

    new_node.type_path_[depth_] = type;
    new_node.deltas_path_[type] = deltas;
    return new_node;
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
