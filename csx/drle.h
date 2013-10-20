/* -*- C++ -*-
 *
 * drle.h -- Delta Run-Length Encoding Manager
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2013, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_DRLE_H__
#define CSX_DRLE_H__

#include "encodings.h"
#include "runtime.h"
#include "sparse_util.h"
#include "sparse_partition.h"
#include "TimingFramework.h"
#include "node.h"
#include "logger.hpp"

#include <map>
#include <bitset>
#include <limits>
#include <cassert>
#include <iostream>
#include <sstream>

#define FOREACH BOOST_FOREACH

using namespace std;

namespace csx {

/**
 *  Delta Run-Length Encoding Manager.
 *
 *  This class is responsible for searching for patterns inside the matrix and
 *  encoding the matrix.
 */
template<typename IndexType, typename ValueType>
class EncodingManager
{
public:    
    /**
     *  Algorithm for spliting the input matrix into sorting windows.
     */
    typedef enum {
        SPLIT_BY_ROWS = 0,  /* Split by rows (each window has the same number of
                               rows) */
        SPLIT_BY_NNZ,       /* Split by nonzeros (each window has almost the
                               same number of nonzero elements) */
    } split_alg_t;
    
    typedef typename std::vector<IndexType>::iterator sort_split_iterator;
    typedef std::map<Encoding::Type, DeltaRLE::Stats> StatsMap;
    
    EncodingManager(SparsePartition<IndexType, ValueType> *spm,
                    RuntimeConfiguration &config);
    ~EncodingManager() {
        if (selected_splits_)
            delete selected_splits_;
    }

    /**
     *  Generate pattern statistics for a sub-matrix of SparsePartition.
     *
     *  @param sp  the matrix to generate statistics for.
     *  @param rs  starting row of matrix where to start searching.
     *  @param re  end row of matrix where to end searching.
     */
    DeltaRLE::Stats GenerateStats(SparsePartition<IndexType, ValueType> *sp,
                                  IndexType rs, IndexType re);
    DeltaRLE::Stats GenerateStats(IndexType rs, IndexType re);
    void GenerateDeltaStats(SparsePartition<IndexType, ValueType> *sp,
                            IndexType rs, IndexType re, DeltaRLE::Stats &stats);

    /**
     *  Generate statistics for all available patterns for the matrix owned by
     *  this EncodingManager.
     */
    void GenAllStats();

    /**
     *  Output statistics to a specified output stream.
     *
     *  @param os the output stream where the statistics will be outputed.
     */
    void OutStats(ostringstream &os);

    /**
     *  Instruct EncodingManager to ignore a specific type of pattern. Patterns
     *  of this type will not be considered during generation of statistics.
     *
     *  @param type the pattern type to be ignored.
     *  @see GenAllStats()
     */
    void AddIgnore(Encoding::Type type)
    {
        xforms_ignore_.set(type);
    }

    void AddIgnore(const Encoding &enc)
    {
        vector<Encoding::Type> types;
        enc.GetTypes(types);
        for (size_t i = 0; i < types.size(); i++)
            AddIgnore(types[i]);
    }

    /**
     *  Instruct EncodingManager to ignore all types of patterns.
     */
    void IgnoreAll()
    {
        xforms_ignore_.set();
    }

    void RemoveIgnore(const Encoding::Type type)
    {
        if (!onedim_blocks_ && (type == Encoding::BlockRow1 ||
                                type == Encoding::BlockCol1))
            return;
        
        assert(!Encoding(type).IsGroup());
        xforms_ignore_.reset(type);
    }

    void RemoveIgnore(const Encoding &enc)
    {
        vector<Encoding::Type> types;
        enc.GetTypes(types);
        for (size_t i = 0; i < types.size(); i++)
            RemoveIgnore(types[i]);
    }

    void RemoveIgnore(const EncodingSequence &encseq)
    {
        for (EncodingSequence::const_iterator iter = encseq.Cbegin();
             iter != encseq.Cend(); ++iter) {
            const Encoding &e = iter->first;
            RemoveIgnore(e);
        }
    }

    /**
     *  Instruct EncodingManager not to ignore any type of available patterns.
     *
     *  @see IgnoreAll()
     */
    void RemoveAll()
    {
        xforms_ignore_.reset();
    }

    /**
     *  Choose a pattern type for encoding the matrix. This function chooses the
     *  type with the highest score.
     *
     *  @return the type chosen.
     *  @see GetTypeScore()
     */
    Encoding::Type ChooseType();

    /**
     *  Get the score of the specified type of patterns. This function
     *  implements a heuristic to characterize the candidate encoding patterns.
     *  Statistics for the different patterns must have already been gathered.
     *  Currently, the heuristic implemented by this function just computes
     *  the following:
     *
     *  <br/><br/>
     *  <i>nr_nzeros_encoded - nr_patterns</i>
     *  <br/><br/>
     *
     *  This is actually a metric of the compression that the specified type of
     *  patterns can achieve.
     *
     *  @param type the type whose score will be evaluated.
     *  @return     score of the type.
     *  @see GenAllStats()
     */
    long GetTypeScore(Encoding::Type type);

    /**
     *  Encode all the patterns of the specified type.
     *
     *  @param type         type of patterns to be encoded.
     */
    void Encode(Encoding::Type type = Encoding::None);

    /**
     *  Decode, i.e., restore to their original form, all the matrix elements of
     *  the specified type.
     *
     *  @param type type of patterns to be decoded.
     */
    void Decode(Encoding::Type type = Encoding::None);

    /**
     *  Encode all available patterns in the matrix owned by this
     *  EncodingManager and send the output (statistics, choice, etc.) to the
     *  specified output stream.
     *
     *  @param os           the output stream, where to send the output.
     *  @see Encode()
     */
    void EncodeAll(ostringstream &os);

    /**
     *  Search for and output all the encoding sequences of the matrix owned
     *  by this EncodingManager.
     */
    void MakeEncodeTree();

    /**
     *  Encode patterns of a series of types explicitly specified by the user.
     */
    void EncodeSerial(const EncodingSequence &encseq);

    /**
     *  Output the split points of the sorting windows, if windows are used for
     *  the encoding of the matrix.
     *
     *  @param out the output stream, where output is to be sent.
     */
    void OutputSortSplits(ostringstream& out);

    void EnableOneDimBlocks(bool enable)
    {
        onedim_blocks_ = enable;
    }

private:
    /**
     *  Encode a sequence of elements. This function calls DoEncodeBlock() or
     *  DoEncodeBlockAlt() if the type of the matrix owned by this
     *  EncodingManager is a block type.
     *
     *  @param xs           the sequence of the column indices of the elements
     *                      to encode. This vector will be cleared at exit.
     *  @param vs           the values of elements to be encoded. This vector
     *                      will be cleared at exit.
     *  @param encoded      vector to append the encoded elements.
     *  @see DoEncodeBlock()
     *  @see DoEncodeBlockAlt()
     */
    void DoEncode(std::vector<IndexType> &xs, std::vector<ValueType> &vs,
                  std::vector<Elem<IndexType, ValueType> > &encoded);

    /**
     *  Encode elements for a block type of patterns.
     *
     *  @param xs      deltas of elements that are going to be encoded.
     *  @param vs      values of elements that are going to be encoded.
     *  @param encoded vector to append the encoded elements.
     */
    void DoEncodeBlock(std::vector<IndexType> &xs, std::vector<ValueType> &vs,
                       std::vector<Elem<IndexType, ValueType> > &encoded);

    /**
     *  Encode elements from a block when the split_blocks_ function is active.
     *
     *  @param xs      deltas of elements that are going to be encoded.
     *  @param vs      values of elements that are going to be encoded.
     *  @param encoded vector to append the encoded elements.
     */
    void DoEncodeBlockAlt(std::vector<IndexType> &xs,
                          std::vector<ValueType> &vs,
                          std::vector<Elem<IndexType, ValueType> > &encoded);

    /**
     *  Decode a CSX element.
     *
     *  @param elem   element to decode.
     *  @param newrow vector to append the decoded elements.
     */
    void DoDecode(const Elem<IndexType, ValueType> *elem,
                  std::vector<Elem<IndexType, ValueType> > &newrow);

    /**
     *  Encode the elements of a row.
     *
     *  @param rstart first element of row.
     *  @param rend   last element of row.
     *  @param newrow vector to append the encoded elements.
     */
    void EncodeRow(const Elem<IndexType, ValueType> *rstart,
                   const Elem<IndexType, ValueType> *rend,
                   std::vector<Elem<IndexType, ValueType> > &newrow);

    /**
     *  Decode the elements of a row.
     *
     *  @param rstart first element of row.
     *  @param rend   last element of row.
     *  @param newrow vector to append the decoded elements.
     */
    void DecodeRow(const Elem<IndexType, ValueType> *rstart,
                   const Elem<IndexType, ValueType> *rend,
                   std::vector<Elem<IndexType, ValueType> > &newrow);

    /**
     *  Cut the blocks that pass the max elements limit.
     *
     *  @param stats       stats of the matrix.
     *  @param block_align block's main dimension.
     */
    void CorrectBlockStats(DeltaRLE::Stats *stats, uint64_t block_align);

    /**
     *  Update stats taking into account the split block function.
     *
     *  @param stats       stats of the matrix.
     *  @param size        number of non-zero elements of initial matrix.
     *  @param block_align block's main dimension.
     */
    void HandleStats(DeltaRLE::Stats *stats, uint64_t size,
                     uint64_t block_align);

    /**
     *  Update the statistics of a sub-matrix.
     *
     *  @param sp    object of SparsePartition class that describes
     *               the sub-matrix elements.
     *  @param xs    deltas of elements examined.
     *  @param stats stats of the sub-matrix so far.
     */
    void UpdateStats(SparsePartition<IndexType, ValueType> *sp,
                     std::vector<IndexType> &xs,
                     DeltaRLE::Stats &stats,
                     std::vector<Elem<IndexType, ValueType> *> &elems);

    /**
     *  Add stats of a type.
     *
     *  @param stats stats added.
     *  @param type  type which stats are being updated.
     */
    void UpdateStats(Encoding::Type type, DeltaRLE::Stats stats);

    /**
     *  Updates the stats for a block pattern.
     *
     *  @param xs    deltas of elements examined.
     *  @param stats stats of the sub-matrix so far.
     *  @param align block's main dimension.
     */
    void UpdateStatsBlock(std::vector<IndexType> &xs,
                          DeltaRLE::Stats &stats, uint64_t align,
                          std::vector<Elem<IndexType, ValueType> *> &elems);

    /**
     *  Used when windows mode is enabled to adapt the stats.
     *
     *  @param type   type which stats are being updated.
     *  @param factor weight of the specific type.
     */
    void CorrectStats(Encoding::Type type, double factor);

    /**
     *  Computes the actual splits points (i.e., the starting rows) of the
     *  sorting windows according to the policy of this EncodingManager.
     *
     *  If the split policy is SPLIT_BY_ROWS, the window size specifies the
     *  number of rows of the original matrix assigned to each window.
     *
     *  If the split policy is SPLIT_BY_NNZ, the window size specifies the
     *  number of nonzeros of the original matrix assigned to each window.  If
     *  the SPLIT_BY_NNZ policy is set, the window size is "rounded" to the next
     *  or previous full row, whichever is closer in terms of nonzero elements.
     */
    void ComputeSortSplits();

    void SelectSplits();

    /*
     *  The actual splitting methods.
     */
    void DoComputeSortSplitsByRows();
    void DoComputeSortSplitsByNNZ();

    /*
     *  Check pre-conditions.
     */
    void CheckAndSetSorting();
    void DoCheckSortByRows();
    void DoCheckSortByNNZ();
    void CheckPropability(double probability)
    {
        assert((probability >= 0.0 && probability <= 1.0) &&
               "invalid sampling probability");
    }

private:
    SparsePartition<IndexType, ValueType> *spm_;
    RuntimeConfiguration &config_;
    long min_limit_;
    long max_limit_;
    double min_perc_;
    bool sampling_enabled_;
    uint64_t sort_window_size_;
    PreprocessingMethod preproc_method_;
    double sampling_portion_;
    uint64_t samples_max_;
    bool symmetric_;
    bool split_blocks_;
    bool onedim_blocks_;
    split_alg_t split_type_;
    std::vector<IndexType> sort_splits_;
    std::vector<IndexType> sort_splits_nzeros_;
    size_t *selected_splits_;
    StatsMap stats_;
    std::map<Encoding::Type, std::set<uint64_t> > deltas_to_encode_;
    std::bitset<Encoding::Max> xforms_ignore_;
    TimingFramework pre_timers_;
};

/**
 *  Output the pattern statistics of a matrix SparsePartition,
 *  specified in stats.
 *
 *  @param stats statistics to output.
 *  @param SparsePartition   the matrix whose statistics will be outputted
 *  @param os    output stream, where the statistics will be sent.
 */
template<typename IndexType, typename ValueType>
void DRLE_OutStats(DeltaRLE::Stats &stats,
                   SparsePartition<IndexType, ValueType> &sp, 
                   ostringstream &os);

}

/* Helper functions */
template<typename T> 
struct RLE {
    long freq;
    T val;
};

//
//  FIXME: MaxDelta() could be possibly merged with DeltaEncode()
// 
template<typename IndexType>
IndexType MaxDelta(std::vector<IndexType> xs)
{
    IndexType prev_x, curr_x, max_delta;

    max_delta = 0;
    prev_x = curr_x = xs[0];
    for (size_t i = 1; i < xs.size(); ++i) {
        curr_x = xs[i];
        IndexType delta = curr_x - prev_x;
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

/*
 * Implementation of class EncodingManager
 */
template<typename IndexType, typename ValueType>
EncodingManager<IndexType, ValueType>::
EncodingManager(SparsePartition<IndexType, ValueType> *spm,
                RuntimeConfiguration &config)
        : spm_(spm),
          config_(config),
          min_limit_(config_.GetProperty<size_t>(
                         RuntimeConfiguration::MatrixMinUnitSize)),
          max_limit_(config_.GetProperty<size_t>(
                         RuntimeConfiguration::MatrixMaxUnitSize)),
          min_perc_(config_.GetProperty<double>(
                        RuntimeConfiguration::MatrixMinCoverage)),
          sort_window_size_(config_.GetProperty<size_t>(
                                RuntimeConfiguration::PreprocWindowSize)),
          preproc_method_(config_.GetProperty<PreprocessingMethod>(
                              RuntimeConfiguration::PreprocMethod)),
          sampling_portion_(config_.GetProperty<double>(
                                RuntimeConfiguration::PreprocSamplingPortion)),
          samples_max_(config_.GetProperty<size_t>(
                           RuntimeConfiguration::PreprocNrSamples)),
          split_blocks_(config_.GetProperty<bool>(
                            RuntimeConfiguration::MatrixSplitBlocks)),
          onedim_blocks_(config_.GetProperty<bool>(
                             RuntimeConfiguration::MatrixOneDimBlocks)),
          split_type_(SPLIT_BY_NNZ)
{
    // Ignore all encodings by default
    IgnoreAll();

    // if (!onedim_blocks_) {
    //     AddIgnore(Encoding::BlockRow1);
    //     AddIgnore(Encoding::BlockCol1);
    // }
    
    // Setup preprocessing
    if (preproc_method_.GetType() == PreprocessingMethod::None) {
        sampling_enabled_ = false;
        selected_splits_ = 0;
    } else {
        assert(samples_max_ > 0 && "invalid samples number");
        assert(sampling_portion_ > 0 && sampling_portion_ <= 1 &&
               "invalid sampling portion");
        if (preproc_method_.GetType() == PreprocessingMethod::FixedPortion) {
            sort_window_size_ =
                sampling_portion_*spm_->GetNrNonzeros() / samples_max_;
        }

        ComputeSortSplits();
        if (samples_max_ > sort_splits_.size())
            samples_max_ = sort_splits_.size();

        SelectSplits();
    }
}
    
template<typename IndexType, typename ValueType>
DeltaRLE::Stats EncodingManager<IndexType, ValueType>::GenerateStats(IndexType rs,
                                                                  IndexType re)
{
    return GenerateStats(spm_, rs, re);
}

template<typename IndexType, typename ValueType>
DeltaRLE::Stats EncodingManager<IndexType, ValueType>::
GenerateStats(SparsePartition<IndexType, ValueType> *sp, IndexType rs,
              IndexType re)
{
    std::vector<IndexType> xs;
    DeltaRLE::Stats stats;
    std::vector<Elem<IndexType, ValueType> *> elems;

    for (IndexType i = rs; i < re; ++i) {
        for (Elem<IndexType, ValueType> *elem = sp->RowBegin(i);
             elem != sp->RowEnd(i); ++elem) {
            if (elem->pattern == NULL) {
                xs.push_back(elem->col);
                elems.push_back(elem);
                continue;
            }

            UpdateStats(sp, xs, stats, elems);
        }

        UpdateStats(sp, xs, stats, elems);
    }

    return stats;
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
GenerateDeltaStats(SparsePartition<IndexType, ValueType> *sp, IndexType rs,
                   IndexType re, DeltaRLE::Stats &stats)
{
    std::vector<IndexType> xs;
    
    for (IndexType i = rs; i < re; ++i) {
        int delta_bytes = 0;
        for (Elem<IndexType, ValueType> *elem = sp->RowBegin(i);
             elem != sp->RowEnd(i); ++elem) {
            if (!elem->in_pattern)
                xs.push_back(elem->col);
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::GenAllStats()
{
    DeltaRLE::Stats::iterator iter, tmp;
    DeltaRLE::Stats *sp;
    
    stats_.clear();
    deltas_to_encode_.clear();

    if (sampling_enabled_ && spm_->GetNrRows() > samples_max_) {
        uint64_t samples_nnz = 0;

        spm_->Transform(Encoding::Horizontal);
        for (size_t i = 0; i < samples_max_; ++i) {
            // uint64_t window_start = sort_splits_[selected_splits_[i]];
            // uint64_t window_size =
            //              sort_splits_[selected_splits_[i]+1] - window_start;
            size_t window_start = sort_splits_[selected_splits_[i]];
            size_t window_size =
                         sort_splits_[selected_splits_[i]+1] - window_start;

            SparsePartition<IndexType, ValueType> *window =
                spm_->GetWindow(window_start, window_size);

            if (window->GetNrNonzeros() != 0) {
                samples_nnz += sort_splits_nzeros_[selected_splits_[i]];

                for (Encoding::Type t = Encoding::Horizontal;
                     t < Encoding::Max; ++t) {
                    if (xforms_ignore_[t])
                        continue;

                    DeltaRLE::Stats w_stats;
                    window->Transform(t);
                    w_stats = GenerateStats(window, 0, window->GetRowptrSize() - 1);
                    UpdateStats(t, w_stats);
                }
            }
            window->Transform(Encoding::Horizontal);
            spm_->PutWindow(window);
            delete window;
        }

        for (Encoding::Type t = Encoding::Horizontal; t < Encoding::Max; ++t) {
            if (xforms_ignore_[t])
                continue;

            Encoding e(t);
            size_t block_align = e.GetBlockAlignment();

            CorrectStats(t, spm_->GetNrNonzeros() / (double) samples_nnz);
            sp = &stats_[t];
            if (block_align && split_blocks_)
                CorrectBlockStats(sp, block_align);

            for (iter = sp->begin(); iter != sp->end(); ) {
                tmp = iter++;
                double p = (double) tmp->second.nnz / (double) spm_->GetNrNonzeros();
                if (block_align == 1 && (unsigned int) tmp->first < min_limit_)
                    // The dimension of one-dimensional blocks must exceed
                    // min_limit_.
                    sp->erase(tmp);
                else if (p < min_perc_ || tmp->first >= CSX_PID_OFFSET)
                    sp->erase(tmp);
                else
                    deltas_to_encode_[t].insert(tmp->first);
            }
        }
    } else {
    // first generate stats for the delta type
#ifdef SPM_HEUR_NEW
        GenerateDeltaStats(spm_, 0, spm_->GetRowptrSize() - 1, stats_[NONE]);
#endif
        for (Encoding::Type t = Encoding::Horizontal; t < Encoding::Max; ++t) {
            if (xforms_ignore_[t])
                continue;

            Encoding e(t);
            size_t block_align = e.GetBlockAlignment();

            spm_->Transform(t);
            stats_[t] = GenerateStats(0, spm_->GetRowptrSize() - 1);
            sp = &stats_[t];
            if (block_align && split_blocks_)
                CorrectBlockStats(sp, block_align);

            for (iter = sp->begin(); iter != sp->end(); ) {
                tmp = iter++;
                double p = (double) tmp->second.nnz / (double) spm_->GetNrNonzeros();
                if (block_align == 1 && (unsigned int) tmp->first < min_limit_)
                    // The dimension of one-dimensional blocks must exceed
                    // min_limit_.
                    sp->erase(tmp);
                else if (p < min_perc_ || tmp->first >= CSX_PID_OFFSET)
                    sp->erase(tmp);
                else
                    deltas_to_encode_[t].insert(tmp->first);
            }
#ifdef SPM_HEUR_NEW
            GenerateDeltaStats(spm_, 0, spm_->GetRowptrSize() - 1, *sp);
#endif
        }
    }
}

namespace csx {

template<typename IndexType, typename ValueType>
void DRLE_OutStats(DeltaRLE::Stats &stats,
                   SparsePartition<IndexType, ValueType> &sp, ostringstream &os)
{
    DeltaRLE::Stats::iterator iter;

    for (iter = stats.begin(); iter != stats.end(); ++iter) {
        os << "    " << iter->first << "-> " << "np:"
           << iter->second.npatterns
           << " nnz: " << iter->second.nnz
           << " ("
           <<  100 * ((double) iter->second.nnz / (double) sp.GetNrNonzeros())
           << "%" << ")";
    }
}

} // end of csx namespace

template<typename IndexType, typename ValueType>
Encoding::Type EncodingManager<IndexType, ValueType>::ChooseType()
{
    Encoding::Type ret;
    long max_score;
    EncodingManager::StatsMap::iterator iter;

    ret = Encoding::None;
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

template<typename IndexType, typename ValueType>
long EncodingManager<IndexType, ValueType>::GetTypeScore(Encoding::Type type)
{
    DeltaRLE::Stats *sp;
    DeltaRLE::Stats::iterator iter;
    long ret;

    if (stats_.find(type) == stats_.end())
        return 0;

    long nr_nzeros_encoded = 0;
    long nr_patterns = 0;
    long nr_deltas = spm_->GetNrDeltas();
    sp = &stats_[type];

    for (iter = sp->begin(); iter != sp->end(); ++iter) {
        if (iter->first == 0 && type != Encoding::None) {
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::Encode(Encoding::Type type)
{
    typename SparsePartition<IndexType, ValueType>::Builder *SpmBld = 0;
    std::vector<Elem<IndexType, ValueType> > new_row;
    uint64_t nr_size;
    Elem<IndexType, ValueType> *elems;

    if (type == Encoding::None && ((type = ChooseType()) == Encoding::None))
        return;

    // Transform matrix to the desired iteration order
    spm_->Transform(type);
    
    SpmBld = new typename SparsePartition<IndexType, ValueType>::
        Builder(spm_,  spm_->GetRowptrSize(), spm_->GetElemsSize());

    for (size_t i = 0; i < spm_->GetRowptrSize() - 1; ++i) {
        EncodeRow(spm_->RowBegin(i), spm_->RowEnd(i), new_row);
        nr_size = new_row.size();
        if (nr_size > 0) {
            //pre_timers_[PRE_TIMER_ALLOC].Start();
            if (pre_timers_.Enabled(PRE_TIMER_ALLOC))
                pre_timers_.timers_[PRE_TIMER_ALLOC].Start();
            elems = SpmBld->AllocElems(nr_size);
            //pre_timers_[PRE_TIMER_ALLOC].Pause();
            if (pre_timers_.Enabled(PRE_TIMER_ALLOC))
                pre_timers_.timers_[PRE_TIMER_ALLOC].Pause();
            for (uint64_t i = 0; i < nr_size; ++i)
                MakeRowElem(new_row[i], elems + i);
        }

        new_row.clear();
        SpmBld->NewRow();
    }

    SpmBld->Finalize();
    delete SpmBld;
    AddIgnore(type);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::Decode(Encoding::Type type)
{
    typename SparsePartition<IndexType, ValueType>::Builder *SpmBld = 0;
    std::vector<Elem<IndexType, ValueType> > new_row;
    uint64_t nr_size;
    Elem<IndexType, ValueType> *elems;

    if (type == Encoding::None)
        return;

    // Transform matrix to the desired iteration order
    spm_->Transform(type);

    // Do the decoding
    SpmBld = new typename SparsePartition<IndexType, ValueType>::
        Builder(spm_,  spm_->GetRowptrSize(), spm_->GetElemsSize());

    for (size_t i = 0; i < spm_->GetRowptrSize() - 1; ++i) {
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
    RemoveIgnore(type);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::EncodeAll(std::ostringstream &os)
{
    Encoding::Type type = Encoding::None;
    Encoding::Type enc_seq[Encoding::Max];
    StatsMap::iterator iter; 
    int counter = 0;

    //pre_timers_[PRE_TIMER_TOTAL].Start();
    if (pre_timers_.Enabled(PRE_TIMER_TOTAL))
        pre_timers_.timers_[PRE_TIMER_TOTAL].Start();
    for (;;) {
        //pre_timers_[PRE_TIMER_STATS].Start();
        if (pre_timers_.Enabled(PRE_TIMER_STATS))
            pre_timers_.timers_[PRE_TIMER_STATS].Start();
        GenAllStats();
        //pre_timers_[PRE_TIMER_STATS].Pause();
        if (pre_timers_.Enabled(PRE_TIMER_STATS))
            pre_timers_.timers_[PRE_TIMER_STATS].Pause();
        OutStats(os);
        type = ChooseType();
        if (type == Encoding::None)
            break;

        //pre_timers_[PRE_TIMER_ENCODE].Start();
        if (pre_timers_.Enabled(PRE_TIMER_ENCODE))
            pre_timers_.timers_[PRE_TIMER_ENCODE].Start();
        Encoding e(type);
        os << "Encode to " << e << "\n";
        
        Encode(type);

        // Add the deltas produced by this encoding
        //spm_->nr_deltas_ += stats_[type][0].npatterns;
        spm_->AddNrDeltas(stats_[type][0].npatterns);

        enc_seq[counter++] = type;
        //pre_timers_[PRE_TIMER_ENCODE].Pause();
        if (pre_timers_.Enabled(PRE_TIMER_ENCODE))
            pre_timers_.timers_[PRE_TIMER_ENCODE].Pause();
    }

    spm_->Transform(Encoding::Horizontal);
    //pre_timers_[PRE_TIMER_TOTAL].Pause();
    if (pre_timers_.Enabled(PRE_TIMER_TOTAL))
        pre_timers_.timers_[PRE_TIMER_TOTAL].Pause();

    os << "Encoding sequence: ";
    if (counter == 0) {
        Encoding e(Encoding::None);
        os << e;
    } else {
        Encoding e(enc_seq[0]);
        os << e;
    }

    for (int i = 1; i < counter; i++) {
        Encoding e(enc_seq[i]);
        os << ", " << e;
    }
    os << "\n";

    pre_timers_.PrintAllStats(os);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::MakeEncodeTree()
{
    uint32_t i,j;
    uint32_t count = 0;
    StatsMap::iterator iter;
    std::vector<Node> nodes;
    std::vector<Node>::iterator it,it2;

    nodes.push_back(Node(0));
    while (nodes.size() != 0) {
        for (i = 0; nodes[0].type_ignore_[i] != Encoding::None; ++i)
            AddIgnore(nodes[0].type_ignore_[i]);
            
        GenAllStats();
        for (i = 0; nodes[0].type_ignore_[i] != Encoding::None; ++i)
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

    LOG_DEBUG << "Tree has " << count << " possible paths.\n";
    exit(0);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::EncodeSerial(
    const EncodingSequence &encseq)
{
    IgnoreAll();
    for (EncodingSequence::const_iterator iter = encseq.Cbegin();
         iter != encseq.Cend(); ++iter) {
        const Encoding &e = iter->first;
        RemoveIgnore(e);
        for (EncodingSequence::delta_const_iterator
                 diter = encseq.DeltasCbegin(iter);
             diter != encseq.DeltasCend(iter); ++diter) {
            deltas_to_encode_[e.GetType()].insert(*diter);
        }

        Encode(e.GetType());
        AddIgnore(e.GetType());
    }

    spm_->Transform(Encoding::Horizontal);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::OutputSortSplits(
    std::ostringstream& out)
{
    sort_split_iterator iter;
    for (iter = sort_splits_.begin(); iter != sort_splits_.end() - 1; ++iter) {
        IndexType rs = *iter;
        IndexType re = *(iter + 1);
        IndexType nnz = spm_->GetRow(re) - spm_->GetRow(rs);

        out << "(rs, re, nnz) = (" << rs << ", " << re << ", " << nnz << ")\n";
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
DoEncode(std::vector<IndexType> &xs, std::vector<ValueType> &vs,
         std::vector<Elem<IndexType, ValueType> > &encoded)
{
    IndexType col;
    std::vector<RLE<IndexType> > rles;
    std::set<uint64_t> *deltas_set;
    typename std::vector<ValueType>::iterator vi = vs.begin();
    Elem<IndexType, ValueType> elem;

    Encoding e(spm_->GetType());
    if (e.IsBlock()) {
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
    deltas_set = &deltas_to_encode_[e.GetType()];

    col = 0;
    elem.pattern = NULL;
    FOREACH (RLE<IndexType> rle, rles) {
        long rle_freq, rle_start;
        rle_freq = rle.freq;
        if (rle_freq != 1 && deltas_set->find(rle.val) != deltas_set->end()) {

            col += rle.val;
            if (col != rle.val) {
                Elem<IndexType, ValueType> *last_elem = &encoded.back();
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
                Elem<IndexType, ValueType> *last_elem;
                long curr_freq = std::min(max_limit_, rle_freq);
                
                elem.col = rle_start;
                encoded.push_back(elem);
                last_elem = &encoded.back();
                last_elem->pattern = new DeltaRLE(curr_freq, rle.val,
                                                  spm_->GetType());
                last_elem->vals = new ValueType[curr_freq];
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
            elem.col = col;
            elem.val = *vi++;
            encoded.push_back(elem);
        }
    }

    assert(vi == vs.end() && "not all elements processed or out of bounds");
    xs.clear();
    vs.clear();
}

//static uint32_t nr_lines = 0;
template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
DoEncodeBlock(std::vector<IndexType> &xs, std::vector<ValueType> &vs,
              std::vector<Elem<IndexType, ValueType> > &encoded)
{
    IndexType col;
    std::vector< RLE<IndexType> > rles;
    std::set<uint64_t> *deltas_set;
    typename std::vector<ValueType>::iterator vi = vs.begin();
    Elem<IndexType, ValueType> elem;
    Elem<IndexType, ValueType> *last_elem;

    // do a delta run-length encoding of the x values
    rles = RLEncode(DeltaEncode(xs));

    // Not all delta rles are to be encoded, only those
    // that are in the ->deltas_to_encode set
    Encoding e(spm_->GetType());
    deltas_set = &deltas_to_encode_[e.GetType()];

    int block_align = e.GetBlockAlignment();
    assert(block_align && "not a block type");

    col = 0;
    elem.pattern = NULL;
    FOREACH (RLE<IndexType> rle, rles) {
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
                elem.col = rle_start + i;
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
                elem.col = rle_start + skip_front + i * nr_elem_block;
                encoded.push_back(elem);
                last_elem = &encoded.back();
                Encoding e(spm_->GetType());
                last_elem->pattern = new BlockRLE(nr_elem_block,
                                                  nr_elem_block / block_align,
                                                  e);
                last_elem->vals = new ValueType[nr_elem_block];
                std::copy(vi, vi + nr_elem_block, last_elem->vals);
                vi += nr_elem_block;
            }

            // Add the remaining elements
            for (uint64_t i = 0; i < skip_back; ++i) {
                elem.col = rle_start + skip_front + nr_elem_block * nr_blocks + i;
                elem.val = *vi++;
                encoded.push_back(elem);
            }
        } else {
            // add individual elements
            for (int i = 0; i < rle.freq; ++i) {
                elem.col = col + i * rle.val;
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
DoEncodeBlockAlt(std::vector<IndexType> &xs, std::vector<ValueType> &vs,
                 std::vector<Elem<IndexType, ValueType> > &encoded)
{
    IndexType col;
    std::vector< RLE<IndexType> > rles;
    std::set<uint64_t> *deltas_set;
    typename std::vector<ValueType>::iterator vi = vs.begin();
    Elem<IndexType, ValueType> elem;
    Elem<IndexType, ValueType> *last_elem;

    // do a delta run-length encoding of the x values
    rles = RLEncode(DeltaEncode(xs));

    // Not all delta rles are to be encoded, only those
    // that are in the ->deltas_to_encode set
    Encoding e(spm_->GetType());

    deltas_set = &deltas_to_encode_[e.GetType()];
    int block_align = e.GetBlockAlignment();
    assert(block_align && "not a block type");

    col = 0;
    elem.pattern = NULL;
    FOREACH (RLE<IndexType> rle, rles) {
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
                elem.col = rle_start++;
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
                    
                    elem.col = rle_start;
                    rle_start += nr_elem_block;
                    encoded.push_back(elem);
                    last_elem = &encoded.back();

                    Encoding be(spm_->GetType());
                    last_elem->pattern = new BlockRLE(nr_elem_block, *i, be);
                    last_elem->vals = new ValueType[nr_elem_block];
                    std::copy(vi, vi + nr_elem_block, last_elem->vals);
                    vi += nr_elem_block;
                    nr_elem -= nr_elem_block;
                    other_dim -= (*i);
                }
            }

            // Add the remaining elements
            skip_back += nr_elem;
            for (uint64_t i = 0; i < skip_back; ++i) {
                elem.col = rle_start++;
                elem.val = *vi++;
                encoded.push_back(elem);
            }
        } else {
            // add individual elements
            for (int i = 0; i < rle.freq; ++i) {
                elem.col = col + i * rle.val;
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
DoDecode(const Elem<IndexType, ValueType> *elem, 
         std::vector<Elem<IndexType, ValueType> > &newrow)
{
    long int i;
    IndexType cur_x;
    Elem<IndexType, ValueType> new_elem;

    new_elem.pattern = NULL;
    cur_x = elem->col;
    for (i = 0; i < elem->pattern->GetSize(); ++i) {
        new_elem.col = cur_x;
        new_elem.val = elem->vals[i];
        newrow.push_back(new_elem);
        cur_x = elem->pattern->GetNextCol(cur_x);
    }

    delete elem->pattern;
    delete elem->vals;  // I beleieve this should be delete[]
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
EncodeRow(const Elem<IndexType, ValueType> *rstart,
          const Elem<IndexType, ValueType> *rend,
          std::vector<Elem<IndexType, ValueType> > &newrow)
{
    std::vector<IndexType> xs;
    std::vector<ValueType> vs;

    // gather x values into xs vector until a pattern is found
    // and encode them using DoEncode()
    for (const Elem<IndexType, ValueType> *e = rstart; e < rend; ++e) {
        if (e->pattern == NULL) {
            xs.push_back(e->col);
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
DecodeRow(const Elem<IndexType, ValueType> *rstart,
          const Elem<IndexType, ValueType> *rend,
          std::vector<Elem<IndexType, ValueType> > &newrow)
{
    for (const Elem<IndexType, ValueType> *e = rstart; e < rend; ++e) {
        if (e->pattern != NULL && e->pattern->GetType() == spm_->GetType()) {
            DoDecode(e, newrow);
        } else {
            newrow.push_back(*e);
            if (e->pattern)
                delete e->pattern;
        }
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
CorrectBlockStats(DeltaRLE::Stats *sp, uint64_t block_align)
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
        if (((double) iter->second.nnz) / ((double) spm_->GetNrNonzeros()) >= 
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
HandleStats(DeltaRLE::Stats *sp, uint64_t size, uint64_t block_align)
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
UpdateStats(Encoding::Type type, DeltaRLE::Stats stats)
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
UpdateStats(SparsePartition<IndexType, ValueType> *spm,
            std::vector<IndexType> &xs,
            DeltaRLE::Stats &stats,
            std::vector<Elem<IndexType, ValueType> *> &elems)
{
    std::vector< RLE<IndexType> > rles;
    Encoding e(spm->GetType());
    size_t block_align = e.GetBlockAlignment();

    if (block_align) {
        EncodingManager::UpdateStatsBlock(xs, stats, block_align, elems);
        return;
    }

    assert(xs.size() == elems.size());
    if (xs.size() == 0)
        return;

    rles = RLEncode(DeltaEncode(xs));
    IndexType col = 0;
    bool last_rle_patt = false; // turn on for encoded units
    typename std::vector<Elem<IndexType, ValueType>*>::iterator ei =
        elems.begin();
    typename std::vector<Elem<IndexType, ValueType>*>::iterator ee = elems.end();

    FOREACH(RLE<IndexType> &rle, rles) {
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
UpdateStatsBlock(std::vector<IndexType> &xs, DeltaRLE::Stats &stats,
                 uint64_t block_align,
                 std::vector<Elem<IndexType, ValueType>*> &elems)
{
    std::vector<RLE<IndexType> > rles;

    assert(block_align && "not a block type");
    if (xs.size() == 0)
        return;
        
    rles = RLEncode(DeltaEncode(xs));

    IndexType unit_start = 0;
    typename std::vector<Elem<IndexType, ValueType>*>::iterator ei =
        elems.begin();
    typename std::vector<Elem<IndexType, ValueType>*>::iterator ee =
        elems.end();

    FOREACH (RLE<IndexType> &rle, rles) {
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
                for (long i = 0; i < nr_elem; ++i) {
                    assert(ei < ee);
                    if (((uint64_t) i) < stats[other_dim].nnz)
//                    if (i < stats[other_dim].nnz)
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
CorrectStats(Encoding::Type type, double factor)
{
    if (stats_.find(type) == stats_.end())
        return;

    DeltaRLE::Stats *l_stats = &stats_[type];
    DeltaRLE::Stats::iterator iter;

    for (iter = l_stats->begin(); iter != l_stats->end(); ++iter) {
        iter->second.nnz = (uint64_t) (iter->second.nnz * factor);
        iter->second.npatterns = (long) (iter->second.npatterns * factor);
        assert(iter->second.nnz <= (uint64_t) spm_->GetNrNonzeros() &&
               "nonzeros of pattern exceed nonzeros of matrix");
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::SelectSplits()
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

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::ComputeSortSplits()
{
    switch (split_type_) {
    case EncodingManager::SPLIT_BY_ROWS:
        DoComputeSortSplitsByRows();
        break;
    case EncodingManager::SPLIT_BY_NNZ:
        DoComputeSortSplitsByNNZ();
        break;
    default:
        assert(false && "unknown split algorithm");
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::CheckAndSetSorting()
{
    switch (split_type_) {
    case EncodingManager::SPLIT_BY_ROWS:
        DoCheckSortByRows();
        break;
    case EncodingManager::SPLIT_BY_NNZ:
        DoCheckSortByNNZ();
        break;
    default:
        assert(false && "unknown split algorithm");
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::DoComputeSortSplitsByRows()
{
    size_t nr_rows = spm_->GetRowptrSize() - 1;
    size_t i;

    for (i = 0; i <= nr_rows; i += sort_window_size_)
        sort_splits_.push_back(i);
        
    if (i > nr_rows && i - nr_rows < sort_window_size_ / 2) {
        sort_splits_.push_back(nr_rows);
    } else {
        sort_splits_.pop_back();
        sort_splits_.push_back(nr_rows);
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::DoComputeSortSplitsByNNZ()
{
    size_t nzeros_cnt;
    size_t nr_rows = spm_->GetRowptrSize() - 1;

    nzeros_cnt = 0;
    sort_splits_.push_back(0);
    for (size_t i = 0; i < nr_rows; ++i) {
        size_t new_nzeros_cnt = 
            nzeros_cnt + spm_->GetRow(i+1) - spm_->GetRow(i);
        if (new_nzeros_cnt < sort_window_size_) {
            nzeros_cnt = new_nzeros_cnt;
        } else {
            sort_splits_.push_back(i + 1);
            sort_splits_nzeros_.push_back(new_nzeros_cnt);
            nzeros_cnt = 0;
        }
    }

    if (nzeros_cnt) {
        IndexType &last_split_nz = sort_splits_nzeros_.back();
        last_split_nz += nzeros_cnt;
        if (nzeros_cnt > sort_window_size_ / 2) {
            sort_splits_.push_back(nr_rows);
        } else {
            sort_splits_.pop_back();
            sort_splits_.push_back(nr_rows);
            
        }    
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::DoCheckSortByRows()
{
    assert(sort_window_size_ <= (size_t)(spm_->GetRowptrSize() - 1) &&
           "invalid sort window");
    if (sampling_portion_ == 0.0 || sort_window_size_ == (spm_->GetRowptrSize() - 1))
        sampling_enabled_ = false;
    else
        sampling_enabled_ = true;
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::DoCheckSortByNNZ()
{
    assert(sort_window_size_ <= spm_->GetElemsSize() && "invalid sort window");
    if (sampling_portion_ == 0 || sort_window_size_ == spm_->GetElemsSize())
        sampling_enabled_ = false;
    else
        sampling_enabled_ = true;
}


template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::OutStats(std::ostringstream &os)
{
    typename EncodingManager::StatsMap::iterator iter;
    for (iter = stats_.begin(); iter != stats_.end(); ++iter){
#if SPM_HEUR_NEW
        bool valid_stat = (iter->first == 0) || (iter->first != 0 
                                                 && iter->second.size() > 1);
#else
        bool valid_stat = !iter->second.empty();
#endif
        if (valid_stat) {
            Encoding e(iter->first);
            os << e << "\ts:" 
               << GetTypeScore(e.GetType());
            DRLE_OutStats(iter->second, *(spm_), os);
            os << "\n";
        }
    }
}

#endif  // CSX_DRLE_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
