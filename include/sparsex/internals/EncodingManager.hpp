/*
 * \file EncodingManager.hpp
 *
 * \brief Delta Run-Length Encoding Manager
 *
 * Copyright (C) 2009-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2013, Vasileios Karakasis
 * Copyright (C) 2010-2012, Theodoros Gkountouvas
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_ENCODING_MANAGER_HPP
#define SPARSEX_INTERNALS_ENCODING_MANAGER_HPP

#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Delta.hpp>
#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/Encodings.hpp>
#include <sparsex/internals/Runtime.hpp>
#include <sparsex/internals/SparsePartition.hpp>
#include <sparsex/internals/Statistics.hpp>
#include <sparsex/internals/TimerCollection.hpp>
#include <sparsex/internals/logger/Logger.hpp>
#include <bitset>
#include <limits>
#include <cassert>
#include <sstream>
#include <boost/container/vector.hpp>
#include <boost/foreach.hpp>

#define FOREACH BOOST_FOREACH

using namespace std;
using namespace sparsex::runtime;

namespace sparsex {
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
    static void CheckParams(RuntimeConfiguration &config);

    /**
     *  Algorithm for spliting the input matrix into sorting windows.
     */
    typedef enum {
        SPLIT_BY_ROWS = 0,  /* Split by rows (each window has the same number of
                               rows) */
        SPLIT_BY_NNZ        /* Split by nonzeros (each window has almost the
                               same number of nonzero elements) */
    } split_alg_t;
    
    typedef typename std::vector<IndexType>::iterator sort_split_iterator;

    EncodingManager(SparsePartition<IndexType, ValueType> *spm,
                    RuntimeConfiguration &config);
    ~EncodingManager() {
        if (selected_splits_)
            delete[] selected_splits_;
    }

    void DoCreateTimers()
    {
        tc_.CreateTimer("Total", "Total");
        tc_.CreateTimer("Stats", " * Statistics");
        tc_.CreateTimer("Encode", " * Encoding");
        tc_.CreateTimer("Alloc", " * Allocations");
    }

    /**
     *  Generate pattern statistics for a sub-matrix of SparsePartition.
     *
     *  @param sp  the matrix to generate statistics for.
     *  @param rs  starting row of matrix where to start searching.
     *  @param re  end row of matrix where to end searching.
     */
    void GenerateStats(SparsePartition<IndexType, ValueType> *sp,
                       IndexType rs, IndexType re,
                       StatsCollection<StatsData> &stats);
    void GenerateStats(IndexType rs, IndexType re,
                       StatsCollection<StatsData> &stats)
    {
        GenerateStats(spm_, rs, re, stats);
    }

    void GenerateDeltaStats(SparsePartition<IndexType, ValueType> *sp,
                            IndexType rs, IndexType re,
                            StatsCollection<StatsData> &stats);

    /**
     *  Generate statistics for all available patterns for the matrix owned by
     *  this EncodingManager.
     */
    void GenAllStats(StatsCollection<StatsData> &stats);

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
        for (EncodingSequence::const_iterator iter = encseq.cbegin();
             iter != encseq.cend(); ++iter) {
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
    Encoding::Type ChooseType(const StatsCollection<StatsData> &stats);

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
    unsigned long GetTypeScore(Encoding::Type type, const StatsData &sdata);

    /**
     *  Encode all the patterns of the specified type.
     *
     *  @param type         type of patterns to be encoded.
     */
    void Encode(Encoding::Type type = Encoding::None);

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
     *  @param row_no       row number of the encoded elements  
     *  @param xs           the sequence of the column indices of the elements
     *                      to encode. This vector will be cleared at exit.
     *  @param vs           the values of elements to be encoded. This vector
     *                      will be cleared at exit.
     *  @param encoded      vector to append the encoded elements.
     *  @see DoEncodeBlock()
     *  @see DoEncodeBlockAlt()
     */
    void DoEncode(IndexType row_no, vector<IndexType> &xs,
                  vector<ValueType> &vs,
                  boost::container::vector<
                      Element<IndexType, ValueType> > &encoded);

    /**
     *  Encode elements for a block type of patterns.
     *
     *  @param row_no  row number of the encoded elements  
     *  @param xs      deltas of elements that are going to be encoded.
     *  @param vs      values of elements that are going to be encoded.
     *  @param encoded vector to append the encoded elements.
     */
    void DoEncodeBlock(IndexType row_no, vector<IndexType> &xs,
                       vector<ValueType> &vs,
                       boost::container::vector<
                           Element<IndexType, ValueType> > &encoded);

    /**
     *  Encode elements from a block when the split_blocks_ function is active.
     *
     *  @param row_no  row number of the encoded elements  
     *  @param xs      deltas of elements that are going to be encoded.
     *  @param vs      values of elements that are going to be encoded.
     *  @param encoded vector to append the encoded elements.
     */
    void DoEncodeBlockAlt(IndexType row_no, vector<IndexType> &xs,
                          vector<ValueType> &vs,
                          boost::container::vector<
                              Element<IndexType, ValueType> > &encoded);

    /**
     *  Encode the elements of a row.
     *
     *  @param rstart first element of row.
     *  @param rend   last element of row.
     *  @param newrow vector to append the encoded elements.
     */
    void EncodeRow(
        typename SparsePartition<IndexType, ValueType>::iterator &rstart,
        typename SparsePartition<IndexType, ValueType>::iterator &rend,
        boost::container::vector<Element<IndexType, ValueType> > &newrow);

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
                     std::vector<Element<IndexType, ValueType> > &elems,
                     StatsCollection<StatsData> &stats);

    /**
     *  Updates the stats for a block pattern.
     *
     *  @param xs    deltas of elements examined.
     *  @param stats stats of the sub-matrix so far.
     *  @param align block's main dimension.
     */
    void UpdateStatsBlock(Encoding::Type type, std::vector<IndexType> &xs,
                          uint64_t align,
                          std::vector<Element<IndexType, ValueType> > &elems,
                          StatsCollection<StatsData> &stats);

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
    size_t min_limit_;
    size_t max_limit_;
    double min_perc_;
    bool sampling_enabled_;
    bool minimize_cost_;
    size_t sort_window_size_;
    PreprocessingMethod preproc_method_;
    PreprocessingHeuristic preproc_heur_;
    double sampling_portion_;
    size_t samples_max_;
    bool symmetric_;
    bool split_blocks_;
    bool onedim_blocks_;
    split_alg_t split_type_;
    std::vector<size_t> sort_splits_;
    std::vector<IndexType> sort_splits_nzeros_;
    size_t *selected_splits_;
    StatsCollection<StatsData> encoded_stats_; // stats for the *encoded* types
    set<Encoding::Instantiation> encoded_inst_;
    bitset<Encoding::Max> xforms_ignore_;
    timing::TimerCollection tc_;
};

/* Helper functions */
template<typename T> 
struct RLE {
    size_t freq;
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
void EncodingManager<IndexType, ValueType>::
CheckParams(RuntimeConfiguration &config)
{
    PreprocessingMethod method = config.GetProperty<PreprocessingMethod>(
        RuntimeConfiguration::PreprocMethod);

    if (method.GetType() == PreprocessingMethod::FixedPortion) {
        double sampling_portion = config.GetProperty<double>(
            RuntimeConfiguration::PreprocSamplingPortion);
        size_t samples_max = config.GetProperty<size_t>(
            RuntimeConfiguration::PreprocNrSamples);
        if (samples_max <= 0) {
            LOG_ERROR << "invalid number of samples\n";
            exit(1);
        }
        if (sampling_portion <= 0 || sampling_portion > 1) {
            LOG_ERROR << "invalid sampling portion\n";
            exit(1);
        }
        LOG_INFO << "Preprocessing method: fixed portion\n";
        LOG_INFO << "* Sampling portion: " << sampling_portion << "\n";
        LOG_INFO << "* Number of samples: " << samples_max << "\n";
    } else if (method.GetType() == PreprocessingMethod::FixedWindow) {
        size_t sort_window_size = config.GetProperty<size_t>(
            RuntimeConfiguration::PreprocWindowSize);
        size_t samples_max = config.GetProperty<size_t>(
            RuntimeConfiguration::PreprocNrSamples);
        if (samples_max <= 0) {
            LOG_ERROR << "invalid number of samples\n";
            exit(1);
        }
        if (sort_window_size <= 0) {
            LOG_ERROR << "invalid window size\n";
            exit(1);
        }
        LOG_INFO << "Preprocessing method: fixed window\n";
        LOG_INFO << "* Window size: " << sort_window_size << "\n";
        LOG_INFO << "* Number of samples: " << samples_max << "\n";
    }
    PreprocessingHeuristic heur = config.GetProperty<PreprocessingHeuristic>(
        RuntimeConfiguration::PreprocHeuristic);
    LOG_INFO << "Compression heuristic selected: " << heur.GetName() << "\n";
}

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
          preproc_heur_(config_.GetProperty<PreprocessingHeuristic>(
                            RuntimeConfiguration::PreprocHeuristic)),
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
    // Enable timers
    DoCreateTimers();
    // Ignore all encodings by default
    IgnoreAll();

    // Setup preprocessing
    if (preproc_method_.GetType() == PreprocessingMethod::None) {
        sampling_enabled_ = false;
        selected_splits_ = 0;
    } else {
        sampling_enabled_ = true;
        assert(samples_max_ > 0 && "invalid samples number");
        if (preproc_method_.GetType() == PreprocessingMethod::FixedPortion) {
            assert(sampling_portion_ > 0 && sampling_portion_ <= 1 &&
                   "invalid sampling portion");
            sort_window_size_ =
                sampling_portion_*spm_->GetNrNonzeros() / samples_max_;
        } else if (preproc_method_.GetType() ==
                   PreprocessingMethod::FixedWindow) {
            assert(sort_window_size_ > 0 && "invalid window size");
        }

        ComputeSortSplits();
        if (samples_max_ > sort_splits_.size())
            samples_max_ = sort_splits_.size();

        SelectSplits();
    }

    if (preproc_heur_.GetType() == PreprocessingHeuristic::MinCost) {
        minimize_cost_ = true;
    } else {
        minimize_cost_ = false;
    }
}
    
template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
GenerateStats(SparsePartition<IndexType, ValueType> *sp, IndexType rs,
              IndexType re, StatsCollection<StatsData> &stats)
{
    std::vector<IndexType> xs;
    std::vector<Element<IndexType, ValueType> > elems;
    for (IndexType i = rs; i < re; ++i) {
        typename SparsePartition<IndexType, ValueType>::iterator ei =
            sp->begin(i);
        typename SparsePartition<IndexType, ValueType>::iterator ee =
            sp->end(i);
        for (; ei != ee; ++ei) {
            Marker &m = (*ei).GetMarker();
            if (!m.IsMarked(PatternMarker::InPattern)) {
                xs.push_back((*ei).GetCol());
                elems.push_back(*ei);
                continue;
            }

            UpdateStats(sp, xs, elems, stats);
        }

        UpdateStats(sp, xs, elems, stats);
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
GenerateDeltaStats(SparsePartition<IndexType, ValueType> *sp,
                   IndexType rs, IndexType re,
                   StatsCollection<StatsData> &stats)
{
    std::vector<IndexType> xs;
    for (IndexType i = rs; i < re; ++i) {
        typename SparsePartition<IndexType, ValueType>::iterator ei =
            sp->begin(i);
        typename SparsePartition<IndexType, ValueType>::iterator ee =
            sp->end(i);
        for (; ei != ee; ++ei) {
            Marker &em = (*ei).GetMarker();
            if (!em.IsMarked(PatternMarker::InPattern))
                xs.push_back((*ei).GetCol());
            else
                em.Unmark(PatternMarker::InPattern);

            if (em.IsMarked(PatternMarker::PatternStart)) {
                // new pattern in the row creates a new delta unit
                if (xs.size()) {
                    size_t delta_size = GetDeltaSize(MaxDelta(xs));
                    size_t npatt = iceil(xs.size(), max_limit_);
                    size_t nnz = ((sizeof(IndexType) - delta_size)*xs.size()) /
                        sizeof(IndexType);
                    if (sp->GetType() == Encoding::None)
                        stats.AppendStats(make_pair(sp->GetType(), 0),
                                          StatsData(npatt, nnz, npatt));
                    else
                        stats.AppendStats(make_pair(sp->GetType(), 0),
                                          StatsData(0, 0, npatt));
                    xs.clear();
                }

                em.Unmark(PatternMarker::PatternStart);
            }


            assert(!em.IsMarked(PatternMarker::PatternStart));
            assert(!em.IsMarked(PatternMarker::InPattern));
        }
        
        if (xs.size()) {
            size_t delta_size = GetDeltaSize(MaxDelta(xs));
            size_t npatt = iceil(xs.size(), max_limit_);
            size_t nnz = ((sizeof(IndexType) - delta_size)*xs.size()) /
                sizeof(IndexType);
            if (sp->GetType() == Encoding::None)
                stats.AppendStats(make_pair(sp->GetType(), 0),
                                  StatsData(npatt, nnz, npatt));
            else
                stats.AppendStats(make_pair(sp->GetType(), 0),
                                  StatsData(0, 0, npatt));
            xs.clear();
        }
    }
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::GenAllStats(
    StatsCollection<StatsData> &stats)
{
    encoded_inst_.clear();
    if (sampling_enabled_ && spm_->GetRowptrSize()-1 > samples_max_) {
        size_t samples_nnz = 0;
        spm_->Transform(Encoding::Horizontal);
        for (size_t i = 0; i < samples_max_; ++i) {
            size_t window_start = sort_splits_[selected_splits_[i]];
            size_t window_size =
                         sort_splits_[selected_splits_[i]+1] - window_start;

            // Quick fix for windows of size 0
            if (window_start >= sort_splits_[selected_splits_[i]+1] - 1)
                break;
            SparsePartition<IndexType, ValueType> *window =
                spm_->GetWindow(window_start, window_size);

            if (window->GetNrNonzeros() != 0) {
                samples_nnz += sort_splits_nzeros_[selected_splits_[i]];

                for (Encoding::Type t = Encoding::Horizontal;
                     t < Encoding::Max; ++t) {
                    if (xforms_ignore_[t])
                        continue;

                    window->Transform(t);
                    GenerateStats(window, 0, window->GetRowptrSize() - 1,
                                  stats);
                }
            }

            window->Transform(Encoding::Horizontal);
            spm_->PutWindow(window);
            delete window;
        }

        // Scale statistics
        StatsDataScaler scaler(spm_->GetNrNonzeros()/(double) samples_nnz);
        stats.ManipulateStats(scaler);

        if (split_blocks_) {
            // Split blocks
            BlockSplitter bsplitter(max_limit_, spm_->GetNrNonzeros(),
                                    min_perc_);
            stats.ManipulateStats(bsplitter);
        }

        // Filter out low coverage instantiations
        CoverageFilter cfilter(spm_->GetNrNonzeros(), min_perc_,
                               encoded_inst_);
        stats.ManipulateStats(cfilter);

    } else {
        if (minimize_cost_) {
            // first generate stats for the delta type
            GenerateDeltaStats(spm_, 0, spm_->GetRowptrSize() - 1, stats);
        }

        for (Encoding::Type t = Encoding::Horizontal; t < Encoding::Max; ++t) {
            if (xforms_ignore_[t])
                continue;

            Encoding e(t);
            size_t block_align = e.GetBlockAlignment();
            spm_->Transform(t);
            GenerateStats(0, spm_->GetRowptrSize() - 1, stats);
            if (block_align && split_blocks_) {
                // Split blocks
                BlockSplitter bsplitter(max_limit_, spm_->GetNrNonzeros(),
                                        min_perc_);
                stats.ManipulateStats(bsplitter);
            }

            // Filter out low coverage instantiations
            CoverageFilter cfilter(spm_->GetNrNonzeros(), min_perc_,
                                   encoded_inst_);
            stats.ManipulateStats(cfilter);
            // stats_[t] = GenerateStats(0, spm_->GetRowptrSize() - 1);
            // sp = &stats_[t];
            // if (block_align && split_blocks_)
            //     CorrectBlockStats(sp, block_align);

            // for (iter = sp->begin(); iter != sp->end(); ) {
            //     tmp = iter++;
            //     double p = (double) tmp->second.nnz / (double) spm_->GetNrNonzeros();
            //     if (block_align == 1 && (unsigned int) tmp->first < min_limit_)
            //         // The dimension of one-dimensional blocks must exceed
            //         // min_limit_.
            //         sp->erase(tmp);
            //     else if (p < min_perc_ || tmp->first >= CSX_PID_OFFSET)
            //         sp->erase(tmp);
            //     else
            //         deltas_to_encode_[t].insert(tmp->first);
            // }
                if (minimize_cost_)
                    GenerateDeltaStats(spm_, 0, spm_->GetRowptrSize() - 1,
                                       stats);
        }
    }
}

template<typename IndexType, typename ValueType>
Encoding::Type EncodingManager<IndexType, ValueType>::ChooseType(
    const StatsCollection<StatsData> &stats)
{
    Encoding::Type ret = Encoding::None;
    unsigned long max_score = 0;
    StatsCollection<StatsData>::const_iterator si = stats.cbegin();
    StatsCollection<StatsData>::const_iterator se = stats.cend();
    for (; si != se; ++si) {
        unsigned long score = GetTypeScore((*si).first, (*si).second);
        if (score == 0) {
            AddIgnore((*si).first);
        } else if (score > max_score) {
            max_score = score;
            ret = (*si).first;
        }
    }

    return ret;
}

template<typename IndexType, typename ValueType>
unsigned long EncodingManager<IndexType, ValueType>::GetTypeScore(
    Encoding::Type type, const StatsData &sdata)
{
    size_t nr_encoded = sdata.GetEncoded();
    size_t nr_patterns = sdata.GetPatterns();
    size_t score;

    if (minimize_cost_) {
        size_t nr_deltas = encoded_stats_.GetGlobalStats().GetDeltas() +
            sdata.GetDeltas();
        size_t nr_switches;
        if (type == Encoding::None)
            nr_switches = nr_patterns;
        else
            nr_switches = nr_patterns + nr_deltas;
        if (nr_encoded < nr_patterns + nr_switches)
            score = 0;
        else
            score = nr_encoded - nr_patterns - nr_switches;
    } else {
        score = nr_encoded - nr_patterns;
    }

    return score;
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::Encode(Encoding::Type type)
{
    boost::container::vector<Element<IndexType, ValueType> > new_row;

    if (type == Encoding::None)
        return;

    // Transform matrix to the desired iteration order
    spm_->Transform(type);

    typename SparsePartition<IndexType, ValueType>::Builder SpmBld
        (spm_, spm_->GetRowptrSize(), spm_->GetElemsSize());

    for (size_t i = 0; i < spm_->GetRowptrSize() - 1; ++i) {
        typename SparsePartition<IndexType, ValueType>::iterator rbegin =
            spm_->begin(i);
        typename SparsePartition<IndexType, ValueType>::iterator rend =
            spm_->end(i);
        EncodeRow(rbegin, rend, new_row);
        size_t nr_size = new_row.size();
        if (nr_size > 0) {
            tc_.StartTimer("Alloc");
            for (size_t i = 0; i < nr_size; ++i)
                SpmBld.AppendElem(new_row[i]);
            tc_.PauseTimer("Alloc");
        }

        new_row.clear();
        SpmBld.NewRow();
    }

    SpmBld.Finalize();
    AddIgnore(type);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::EncodeAll(ostringstream &os)
{
    if (!spm_->GetNrNonzeros())
        // If matrix is empty, just return
        return;

    Encoding::Type type = Encoding::None;
    Encoding::Type enc_seq[Encoding::Max];
    size_t cnt = 0;
    
    encoded_stats_.Clear();
    tc_.StartTimer("Total");
    for (;;) {
        StatsCollection<StatsData> type_stats;
        tc_.StartTimer("Stats");
        GenAllStats(type_stats);
        tc_.PauseTimer("Stats");
        os << type_stats << "\n";
        type = ChooseType(type_stats);
        if (type == Encoding::None)
            break;

        tc_.StartTimer("Encode");
        Encoding e(type);
        os << "Encode to " << e << "\n";
        
        encoded_stats_.AppendStats(type, type_stats);
        Encode(type);

        enc_seq[cnt++] = type;
        tc_.PauseTimer("Encode");
    }

    spm_->Transform(Encoding::Horizontal);
    tc_.PauseTimer("Total");

    os << "Encoding sequence: ";
    if (!cnt) {
        Encoding e(Encoding::None);
        os << e;
    } else {
        Encoding e(enc_seq[0]);
        os << e;
    }

    for (size_t i = 1; i < cnt; i++) {
        Encoding e(enc_seq[i]);
        os << ", " << e;
    }

    os << "\n";
    os << "==== PREPROCESSING TIMING STATISTICS ====\n";
    tc_.PrintAllTimers(os);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::EncodeSerial(
    const EncodingSequence &encseq)
{
    if (!spm_->GetNrNonzeros())
        // If matrix is empty, just return
        return;

    IgnoreAll();
    for (EncodingSequence::const_iterator iter = encseq.cbegin();
         iter != encseq.cend(); ++iter) {
        const Encoding &e = iter->first;
        RemoveIgnore(e);
        for (EncodingSequence::delta_const_iterator
                 diter = encseq.DeltasCbegin(iter);
             diter != encseq.DeltasCend(iter); ++diter) {
            encoded_inst_.insert(make_pair(e.GetType(), *diter));
        }

        Encode(e.GetType());
        AddIgnore(e.GetType());
    }

    spm_->Transform(Encoding::Horizontal);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::OutputSortSplits(ostringstream& out)
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
DoEncode(IndexType row_no,
         std::vector<IndexType> &xs, std::vector<ValueType> &vs,
         boost::container::vector<Element<IndexType, ValueType> > &encoded)
{
    typename std::vector<ValueType>::iterator vi = vs.begin();

    Encoding e(spm_->GetType());
    if (e.IsBlock()) {
        if (!split_blocks_)
            DoEncodeBlock(row_no, xs, vs, encoded);
        else
            DoEncodeBlockAlt(row_no, xs, vs, encoded);
        return;
    }

    // do a delta run-length encoding of the x values
    std::vector<RLE<IndexType> > rles = RLEncode(DeltaEncode(xs));

    IndexType col = 0;
    FOREACH (RLE<IndexType> rle, rles) {
        size_t rle_freq, rle_start;
        rle_freq = rle.freq;
        if (rle_freq != 1 &&
            encoded_inst_.find(
                make_pair(e.GetType(), rle.val)) != encoded_inst_.end()) {

            col += rle.val;
            if (col != rle.val) {
#if !SPX_USE_NUMA
                Element<IndexType, ValueType> &last_elem = encoded.back();
#endif
                rle_start = col;
                rle_freq = rle.freq;
#if !SPX_USE_NUMA
                if (!last_elem.IsPattern()) {
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
                size_t curr_freq = min(max_limit_, rle_freq);
                encoded.emplace_back(row_no, rle_start, vi, curr_freq,
                                     make_pair(spm_->GetType(), rle.val));

                vi += curr_freq;
                rle_start += rle.val * curr_freq;
                rle_freq -= curr_freq;
                ++nr_parts;
            }

            // leave col at the last element of the pattern
            col = rle_start - rle.val;
        }

        // add individual elements
        for (size_t i = 0; i < rle_freq; ++i) {
            col += rle.val;
            encoded.emplace_back(row_no, col, *vi++);
        }
    }

    assert(vi == vs.end() && "not all elements processed or out of bounds");
    xs.clear();
    vs.clear();
}

//static uint32_t nr_lines = 0;
template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
DoEncodeBlock(IndexType row_no,
              std::vector<IndexType> &xs, std::vector<ValueType> &vs,
              boost::container::vector<Element<IndexType, ValueType> > &encoded)
{
    typename std::vector<ValueType>::iterator vi = vs.begin();

    // do a delta run-length encoding of the x values
    std::vector< RLE<IndexType> > rles = RLEncode(DeltaEncode(xs));

    Encoding e(spm_->GetType());

    int block_align = e.GetBlockAlignment();
    assert(block_align && "not a block type");

    IndexType col = 0;
    FOREACH (RLE<IndexType> rle, rles) {
        size_t skip_front, skip_back, nr_elem;
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
            encoded_inst_.find(
                make_pair(e.GetType(),
                          nr_elem / block_align)) != encoded_inst_.end() &&
            nr_elem >= (size_t) 2 * block_align) {

            size_t rle_start;
            if (col != 1) {
                rle_start = col - 1;
                encoded.pop_back();
                --vi;
            } else {
                rle_start = col;
            }

            // Add elements skipped from start
            for (size_t i = 0; i < skip_front; ++i) {
                encoded.emplace_back(row_no, rle_start + i, *vi++);
            }

            // Align max_limit
            size_t max_limit = max_limit_ / block_align * block_align;
            size_t nr_blocks = nr_elem / max_limit;
            size_t nr_elem_block = min(max_limit, nr_elem);

            if (nr_blocks == 0)
                nr_blocks = 1;
            else
                // Adjust skip_back with leftover items from max_limit alignment
                skip_back += nr_elem - nr_elem_block * nr_blocks;

            for (size_t i = 0; i < nr_blocks; ++i) {
                // Add the blocks
                encoded.emplace_back(row_no,
                                     rle_start + skip_front + i*nr_elem_block,
                                     vi, nr_elem_block,
                                     make_pair(spm_->GetType(),
                                               nr_elem_block / block_align));

                vi += nr_elem_block;
            }

            // Add the remaining elements
            for (size_t i = 0; i < skip_back; ++i) {
                encoded.emplace_back(row_no, rle_start + skip_front +
                                     nr_elem_block * nr_blocks + i,
                                     *vi++);
            }
        } else {
            // add individual elements
            for (size_t i = 0; i < rle.freq; ++i) {
                encoded.emplace_back(row_no, col + i*rle.val, *vi++);
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
DoEncodeBlockAlt(
    IndexType row_no, vector<IndexType> &xs, vector<ValueType> &vs,
    boost::container::vector<Element<IndexType, ValueType> > &encoded)
{
    typename std::vector<ValueType>::iterator vi = vs.begin();

    // do a delta run-length encoding of the x values
    std::vector< RLE<IndexType> > rles = RLEncode(DeltaEncode(xs));

    Encoding e(spm_->GetType());

    size_t block_align = e.GetBlockAlignment();
    assert(block_align && "not a block type");

    IndexType col = 0;
    FOREACH (RLE<IndexType> rle, rles) {
        size_t skip_front, skip_back, nr_elem;

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
        if (rle.val == 1 && nr_elem >= 2*block_align) {
            size_t rle_start;
            if (col != 1) {
                rle_start = col - 1;
                encoded.pop_back();
                --vi;
            } else {
                rle_start = col;
            }

            // Add elements skipped from start
            for (size_t i = 0; i < skip_front; ++i) {
                encoded.emplace_back(row_no, rle_start++, *vi++);
            }

            // Split blocks to match the appropriate sizes for encode
            size_t other_dim = nr_elem / block_align;

            for (set<Encoding::Instantiation>::reverse_iterator i =
                     encoded_inst_.rbegin(); i != encoded_inst_.rend(); ++i) {
                if (i->first != e.GetType())
                    continue;

                while (other_dim >= i->second) {
                    // New block RLE
                    size_t nr_elem_block = block_align * i->second;
                    encoded.emplace_back(row_no, rle_start, vi, nr_elem_block,
                                         make_pair(spm_->GetType(), i->second));

                    rle_start += nr_elem_block;
                    vi += nr_elem_block;
                    nr_elem -= nr_elem_block;
                    other_dim -= i->second;
                }
            }

            // Add the remaining elements
            skip_back += nr_elem;
            for (size_t i = 0; i < skip_back; ++i) {
                encoded.emplace_back(row_no, rle_start++, *vi++);
            }
        } else {
            // add individual elements
            for (size_t i = 0; i < rle.freq; ++i)
                encoded.emplace_back(row_no, col + i*rle.val, *vi++);
        }

        col += rle.val * (rle.freq - 1);
    }

    assert(vi == vs.end() && "out of bounds");
    xs.clear();
    vs.clear();
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
EncodeRow(typename SparsePartition<IndexType, ValueType>::iterator &rstart,
          typename SparsePartition<IndexType, ValueType>::iterator &rend,
          boost::container::vector<Element<IndexType, ValueType> > &newrow)
{
    std::vector<IndexType> xs;
    std::vector<ValueType> vs;

    // gather x values into xs vector until a pattern is found
    // and encode them using DoEncode()
    IndexType row_no = (*rstart).GetRow();
    typename SparsePartition<IndexType, ValueType>::iterator &ei = rstart;
    for (; ei != rend; ++ei) {
        if (!(*ei).IsPattern()) {
            xs.push_back((*ei).GetCol());
            vs.push_back((*ei).GetValue());
            continue;
        }

        if (xs.size() != 0)
            DoEncode(row_no, xs, vs, newrow);

        newrow.push_back(*ei);
    }

    // Encode any remaining elements
    if (xs.size() != 0)
        DoEncode(row_no, xs, vs, newrow);
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
UpdateStats(SparsePartition<IndexType, ValueType> *spm,
            std::vector<IndexType> &xs,
            std::vector<Element<IndexType, ValueType> > &elems,
            StatsCollection<StatsData> &stats)
{
    std::vector< RLE<IndexType> > rles;
    Encoding e(spm->GetType());
    size_t block_align = e.GetBlockAlignment();

    if (block_align) {
        UpdateStatsBlock(spm->GetType(), xs, block_align, elems, stats);
        return;
    }

    assert(xs.size() == elems.size());
    if (xs.size() == 0)
        return;

    rles = RLEncode(DeltaEncode(xs));
    IndexType col = 0;
#if !SPX_USE_NUMA
    bool last_rle_patt = false; // turn on for encoded units
#endif
    typename std::vector<Element<IndexType, ValueType> >::iterator ei =
        elems.begin();

    FOREACH(RLE<IndexType> &rle, rles) {
#if SPX_USE_NUMA
        size_t real_limit = min_limit_;
#else
        size_t real_limit = (col && !last_rle_patt) ?
            min_limit_ - 1 : min_limit_;
#endif
        if (rle.freq > 1 && rle.freq >= real_limit) {
#if SPX_USE_NUMA
            size_t real_nnz = rle.freq;
#else
            size_t real_nnz = (col && !last_rle_patt) ? rle.freq + 1 : rle.freq;
#endif
            size_t rem_nnz = real_nnz % max_limit_;
            size_t patt_nnz = real_nnz;
            size_t patt_npatterns = real_nnz / max_limit_ + (rem_nnz != 0);
            if (rem_nnz && rem_nnz < min_limit_) {
                // remove the pattern for the remainder elements; too short
                --patt_npatterns;
                patt_nnz -= rem_nnz;
            } else {
                rem_nnz = 0;
            }

            stats.AppendStats(make_pair(spm->GetType(), rle.val),
                              StatsData(patt_nnz, patt_npatterns));

            if (real_nnz > rle.freq)
                --ei;

            // Mark the elements of the pattern and its start
            (*ei).GetMarker().Mark(PatternMarker::PatternStart);
            for (size_t i = 0; i < real_nnz - rem_nnz; ++i) {
                assert(ei < elems.end());
                (*ei).GetMarker().Mark(PatternMarker::InPattern);
                ++ei;
            }

#if !SPX_USE_NUMA
            last_rle_patt = true;
#endif
        } else {
#if !SPX_USE_NUMA
            last_rle_patt = false;
#endif
            ei += rle.freq;
        }

        col += rle.val;
    }

    xs.clear();
    elems.clear();
}

template<typename IndexType, typename ValueType>
void EncodingManager<IndexType, ValueType>::
UpdateStatsBlock(Encoding::Type type, std::vector<IndexType> &xs,
                 size_t block_align,
                 std::vector<Element<IndexType, ValueType> > &elems,
                 StatsCollection<StatsData> &stats)
{
    std::vector<RLE<IndexType> > rles;

    assert(block_align && "not a block type");
    if (xs.size() == 0)
        return;
        
    rles = RLEncode(DeltaEncode(xs));

    IndexType unit_start = 0;
    // typename std::vector<Element<IndexType, ValueType>*>::iterator ei =
    //     elems.begin();
    typename std::vector<Element<IndexType, ValueType> >::iterator ei =
        elems.begin();

    FOREACH (RLE<IndexType> &rle, rles) {
        unit_start += rle.val;
        if (rle.val == 1) {
            // Start of the real block is at `unit_start - 1' with
            // one-based indexing. When computing the `%' we need
            // zero-based indexing.
            size_t nr_elem;
            size_t skip_front;

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

            size_t other_dim = nr_elem / block_align;

            ei += min(skip_front, nr_elem);
            if (other_dim >= 2) {
                size_t patt_nnz = other_dim * block_align;
                size_t patt_npatterns = 1;
                stats.AppendStats(make_pair(type, other_dim),
                                  StatsData(patt_nnz, patt_npatterns));

                // Mark the elements
                (*ei).GetMarker().Mark(PatternMarker::PatternStart);
                for (size_t i = 0; i < nr_elem; ++i) {
                    assert(ei < elems.end());
                    if (i < patt_nnz)
                        (*ei).GetMarker().Mark(PatternMarker::InPattern);
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
    if (sampling_portion_ == 0.0 ||
        sort_window_size_ == (spm_->GetRowptrSize() - 1))
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

}   // end of namespace csx
}   // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_ENCODING_MANAGER_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
