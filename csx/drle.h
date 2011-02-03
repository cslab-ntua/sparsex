/* -*- C++ -*-
 *
 * drle.h -- Delta Run-Length Encoding Manager
 *
 * Copyright (C) 2010, Computing Systems Laboratory (CSLab), NTUA.
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_DRLE_H__
#define CSX_DRLE_H__

#include <map>
#include <set>
#include <bitset>
#include <limits>
#include <cassert>
#include <cstdlib>

#include "spm.h"

namespace csx {

/**
 *  Delta Run-Length Encoding Manager.
 *
 *  This class is responsible for searching for patters inside the matrix and
 *  encoding the matrix.
 */
class DRLE_Manager
{
public:
    SPM *spm;
    long min_limit;
    long max_limit;
    double min_perc;

    /**
     *  Algorithm for spliting the input matrix into sorting windows.
     */
    typedef enum {
        SPLIT_BY_ROWS = 0,  ///< Split by rows (each window has the same number
                            ///  of rows)
        SPLIT_BY_NNZ,       ///< Split by nonzeros (each window has almost the
                            ///  same number of nonzero elements)
    } split_alg_t;

    DRLE_Manager(SPM *_spm,
                 long min_limit_ = 4,
                 long max_limit_ = std::numeric_limits<long>::max(),
                 double min_perc_ = .1,
                 uint64_t sort_window_size_ = 0,
                 split_alg_t split_type_ = SPLIT_BY_ROWS,
                 double probability = 1.0,
                 uint64_t samples_max_ = std::numeric_limits<uint64_t>::max());

    /**
     *  Generate pattern statistics for a sub-matrix of spm.
     *
     *  @param spm the matrix to generate statistics for.
     *  @param rs  starting row of matrix where to start searching.
     *  @param re  end row of matrix where to end searching.
     */
    DeltaRLE::Stats GenerateStats(SPM *spm, uint64_t rs, uint64_t re);
    DeltaRLE::Stats GenerateStats(uint64_t rs, uint64_t re);

    typedef std::map<SpmIterOrder, DeltaRLE::Stats> StatsMap;
    StatsMap stats;

    /**
     *  Generate statistics for all available patterns for the matrix owned by
     *  this DRLE_Manager.
     *
     *  @param split_blocks determines whether or not to use the extra function
     *                      of spliting blocks.
     */
    void GenAllStats(bool split_blocks);

    /**
     *  Output statistics to a specified output stream.
     *
     *  @param os the output stream where the statistics will be outputed.
     */
    void OutStats(std::ostream &os=std::cout);

    std::map<SpmIterOrder, std::set<uint64_t> > deltas_to_encode;
    std::bitset<XFORM_MAX> xforms_ignore;

    /**
     *  Instruct DRLE_Manager to ignore a specific type of pattern. Patterns of
     *  this type will not be considered during generation of statistics.
     *
     *  @param type the pattern type to be ignored.
     *  @see GenAllStats()
     */
    void AddIgnore(SpmIterOrder type);

    /**
     *  Instruct DRLE_Manager to ignore all types of patterns.
     */
    void IgnoreAll();

    /**
     *  Instruct DRLE_Manager not to ignore the patterns of type <tt>type</tt>
     *  anymore.
     *
     *  @param type type which is no longer ignored.
     *  @see AddIgnore()
     */
    void RemoveIgnore(SpmIterOrder type);

    /**
     *  Instruct DRLE_Manager not to ignore any type of available patterns.
     *
     *  @see IgnoreAll()
     */
    void RemoveAll();

    /**
     *  Choose a pattern type for encoding the matrix. This function chooses the
     *  type with the highest score.
     *
     *  @see GetTypeScore()
     */
    SpmIterOrder ChooseType();

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
    uint64_t GetTypeScore(SpmIterOrder type);

    /**
     *  Encode all the patterns of the specified type.
     *
     *  @param type         type of patterns to be encoded.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void Encode(SpmIterOrder type = NONE, bool split_blocks = false);

    /**
     *  Decode, i.e., restore to their original form, all the matrix elements of
     *  the specified type.
     *
     *  @param type type of patterns to be decoded.
     */
    void Decode(SpmIterOrder type = NONE);

    /**
     *  Encode all available patterns in the matrix owned by this DRLE_Manager
     *  and send the output (statistics, choice, etc.) to the specified output
     *  stream.
     *
     *  @param os           the output stream, where to send the output.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     *  @see Encode()
     */
    void EncodeAll(std::ostream &os, bool split_blocks);

    /**
     *  Search for and output all the encoding sequences of the matrix owned
     *  by this DRLE_Manager.
     *
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void MakeEncodeTree(bool operate);

    /**
     *  Encode patterns of a series of types explicitly specified by the user.
     *
     *  @param xform        array of types_id.
     *  @param deltas       deltas that match to the corresponding types.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void EncodeSerial(int *xform, int **deltas, bool operate);

    /**
     *  Output the split points of the sorting windows, if windows are used for
     *  the encoding of the matrix.
     *
     *  @param out the output stream, where output is to be sent.
     */
    void OutputSortSplits(std::ostream& out);

private:
    bool sort_windows;
    uint64_t sort_window_size;
    split_alg_t split_type;
    std::vector<uint64_t> sort_splits;
    typedef std::vector<uint64_t>::iterator sort_split_iterator;
    double sampling_probability;
    uint64_t samples_max;
    static const uint64_t max_sampling_tries = 3;

    /**
     *  Encode a sequence of elements. This function calls DoEncodeBlock() or
     *  DoEncodeBlockAlt() if the type of the matrix owned by this DRLE_Manager
     *  is a block type.
     *
     *  @param xs           the sequence of the column indices of the elements
     *                      to encode. This vector will be cleared at exit.
     *  @param vs           the values of elements to be encoded. This vector
     *                      will be cleared at exit.
     *  @param encoded      vector to append the encoded elements.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     *  @see DoEncodeBlock()
     *  @see DoEncodeBlockAlt()
     */
    void DoEncode(std::vector<uint64_t> &xs, std::vector<double> &vs,
                  std::vector<SpmRowElem> &encoded, bool split_blocks);

    /**
     *  Encode elements for a block type of patterns.
     *
     *  @param xs      deltas of elements that are going to be encoded.
     *  @param vs      values of elements that are going to be encoded.
     *  @param encoded vector to append the encoded elements.
     */
    void DoEncodeBlock(std::vector<uint64_t> &xs, std::vector<double> &vs,
                       std::vector<SpmRowElem> &encoded);

    /**
     *  Encode elements from a block when the split_blocks function is active.
     *
     *  @param xs      deltas of elements that are going to be encoded.
     *  @param vs      values of elements that are going to be encoded.
     *  @param encoded vector to append the encoded elements.
     */
    void DoEncodeBlockAlt(std::vector<uint64_t> &xs, std::vector<double> &vs,
                          std::vector<SpmRowElem> &encoded);
    /**
     *  Decode a CSX element.
     *
     *  @param elem   element to decode.
     *  @param newrow vector to append the decoded elements.
     */
    void DoDecode(const SpmRowElem *elem, std::vector<SpmRowElem> &newrow);

    /**
     *  Encode the elements of a row.
     *
     *  @param rstart first element of row.
     *  @param rend   last element of row.
     *  @param newrow vector to append the encoded elements.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void EncodeRow(const SpmRowElem *rstart, const SpmRowElem *rend,
                   std::vector<SpmRowElem> &newrow, bool split_blocks);

    /**
     *  Decode the elements of a row.
     *
     *  @param rstart first element of row.
     *  @param rend   last element of row.
     *  @param newrow vector to append the decoded elements.
     */
    void DecodeRow(const SpmRowElem *rstart, const SpmRowElem *rend,
                   std::vector<SpmRowElem> &newrow);

    /**
     *  Cut the blocks that pass the max elements limit.
     *
     *  @param stats       stats of the matrix.
     *  @param block_align block's main dimension.
     */
    void CutMaxLimit(DeltaRLE::Stats *stats, uint64_t block_align);

    /**
     *  Update stats taking into account the split block function.
     *
     *  @param stats       stats of the matrix.
     *  @param size        number of non-zero elementw of initial matrix.
     *  @param block_align block's main dimension.
     */
    void HandleStats(DeltaRLE::Stats *stats, uint64_t size,
                     uint64_t block_align);

    /**
     *  Update the statistics of a sub-matrix.
     *
     *  @param spm   object of spm class that describes the sub-matrix elements.
     *  @param xs    deltas of elements examined.
     *  @param stats stats of the sub-matrix so far.
     */
    void UpdateStats(SPM *spm, std::vector<uint64_t> &xs,
                     DeltaRLE::Stats &stats);
    void UpdateStats(std::vector<uint64_t> &xs, DeltaRLE::Stats &stats);

    /**
     *  Add stats of a type.
     *
     *  @param stats stats added.
     *  @param type  type which stats are being updated.
     */
    void UpdateStats(SpmIterOrder type, DeltaRLE::Stats stats);

    /**
     *  Updates the stats for a block pattern.
     *
     *  @param xs    deltas of elements examined.
     *  @param stats stats of the sub-matrix so far.
     *  @param align block's main dimension.
     */
    void UpdateStatsBlock(std::vector<uint64_t> &xs,
                          DeltaRLE::Stats &stats, uint64_t align);

    /**
     *  Used when windows mode is enabled to adapt the stats.
     *  @param type   type which stats are being updated.
     *  @param factor weight of the specific type.
     */
    void CorrectStats(SpmIterOrder type, double factor);

    /**
     *  Computes the actual splits points (i.e., the starting rows) of the
     *  sorting windows according to the policy of this DRLE_Manager.
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
};

/**
 *  Output the pattern statistics of a matrix spm, specified in stats.
 *
 *  @param stats statistics to output.
 *  @param spm   the matrix whose statistics will be outputted
 *  @param os    output stream, where the statistics will be sent.
 */
void DRLE_OutStats(DeltaRLE::Stats &stats, SPM &spm, std::ostream &os);

/**
 *  Keeps data of an encoding sequence.
 */
class Node {
public:
    uint32_t depth;
    std::map<SpmIterOrder, std::set<uint64_t> > deltas_path;
    SpmIterOrder *type_path;
    SpmIterOrder *type_ignore;

    Node(uint32_t depth_) : depth(depth_) {
        uint32_t i;
        this->type_path = new SpmIterOrder[depth];
        this->type_ignore = new SpmIterOrder[XFORM_MAX];
        for (i=0; i<((uint32_t) XFORM_MAX); i++)
            this->type_ignore[i] = NONE;
    }
    ~Node() {}
    void PrintNode();

    /**
     *  Ignore the type for the encoding sequence examined.
     *  
     *  @param type type which is ignored.
     */
    void Ignore(SpmIterOrder type);

    /**
     *  Copies a node to a new one and inserts an extra type in the end of the
     *  encoding sequence.
     *  
     *  @param type   type which is inserted in the end.
     *  @param deltas deltas corresponding to type inserted.
     */
    Node MakeChild(SpmIterOrder type, std::set<uint64_t> deltas);
};

}

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
