#ifndef CSX_DRLE_H__
#define CSX_DRLE_H__

#include <map>
#include <set>
#include <bitset>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <stdexcept>

#include "spm.h"

namespace csx {

/**
 *  Responsible for gathering stats, printing stats, encoding and decoding the 
 *  matrix.
 */
class DRLE_Manager
{
public:
    SPM *spm;
    long min_limit;
    long max_limit;
    double min_perc;
    
    typedef enum {
        SPLIT_BY_ROWS = 0,
        SPLIT_BY_NNZ,
    } split_alg_t;
    
    DRLE_Manager(SPM *_spm, long min_limit_=4, 
                 long max_limit_ = std::numeric_limits<long>::max(),
	         double min_perc_=.1, uint64_t sort_window_size_ = 0,
	         split_alg_t split_type_ = SPLIT_BY_ROWS,
	         double probability = 1.0,
	         uint64_t samples_max_ = std::numeric_limits<uint64_t>::max());

    /**
     *  Generate stats for sub-matrix.
     *  @param spm spm class object that describes the sub-matrix.
     *  @param rs  start row of the block in which stats are gathered.
     *  @param re  end row of the block in which stats are gathered.
     */
    DeltaRLE::Stats GenerateStats(SPM *spm, uint64_t rs, uint64_t re);
    DeltaRLE::Stats GenerateStats(uint64_t rs, uint64_t re);
     
    typedef std::map <SpmIterOrder, DeltaRLE::Stats> StatsMap;
    StatsMap stats;
	
    /**
     *  Generates stats for the whole matrix.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void GenAllStats(bool split_blocks);
    
    //TODO : Ask for PrintStats.
    /**
     *  Prints stats.
     *  @param os output channel.
     */
    void OutStats(std::ostream &os=std::cout);
     
    std::map <SpmIterOrder, std::set<uint64_t> > deltas_to_encode;
    std::bitset<XFORM_MAX> xforms_ignore;

    /**
     *  Ignore the specific type. Patterns of this type are no longer available
     *  for gathering stats, encoding or decoding.
     *  @param type type which is ignored.
     */
    void AddIgnore(SpmIterOrder type);
    
    /**
     *  Ignore all the types of patterns.
     */
    void IgnoreAll();
    
    /**
     *  Remove ignore from the specific type (if it has any). Patterns of this
     *  type are now available for gathering stats, encoding or decoding.
     *  @param type type which is no longer ignored.
     */
    void RemoveIgnore(SpmIterOrder type);
    
    /**
     *  Remove ignore for every type of patterns.
     */
    void RemoveAll();

    /**
     *  Heuristic method that chooses the type of patterns that are going to be
     *  encoded. The patterns of the type that is chosen cover larger part of
     *  the matrix than any other type.
     */
    SpmIterOrder ChooseType();
    
    /**
     *  Calculate score of the specific type of patterns. More specifically the
     *  number of elements covered by patterns of this type minus the number of
     *  the patterns is calculated.
     *  @param type type which is going to be evaluated.
     *  @return     score of the specific type.
     */
    uint64_t GetTypeScore(SpmIterOrder type);

    /**
     *  Encode all the patterns of the specific type.
     *  @param type         type that it's patterns are going to be encoded.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void Encode(SpmIterOrder type=NONE, bool split_blocks=false);
	
    /**
     *  Decode all the CSX elements of the specific type.
     *  @param type type that it's CSX elements are going to be decoded.
     */
    void Decode(SpmIterOrder type=NONE);
	
    /**
     *  Gather stats, print them, choose the appropriate type of patterns and 
     *  encode them. Repeat until there is no satisfactory encoding.
     *  @param os           output channel.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void EncodeAll(std::ostream &os, bool split_blocks);
    
    /**
     *  Print all the encoding sequences.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void MakeEncodeTree(bool operate);
    
    /**
     *  Encode patterns of types given by the user.
     *  @param xform        array of types_id.
     *  @param deltas       deltas that match to the corresponding types.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void EncodeSerial(int *xform, int **deltas, bool operate);
    
    void PrintSortSplits(std::ostream& out);

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
     *  Encode some elements from the same row.
     *  @param xs           deltas of elements that are going to be encoded.
     *  @param vs           values of elements that are going to be encoded.
     *  @param newrow       vector to append the encoded elements.
     *  @param split_blocks determines whether or not the extra function of
     *                      spliting blocks is used.
     */
    void DoEncode(std::vector<uint64_t> &xs, std::vector<double> &vs,
                  std::vector<SpmRowElem> &newrow, bool split_blocks);
                      
    //TODO : Merge the two functions DoEncodeBlock into one
    /**
     *  Encode elements from a block.
     *  @param xs     deltas of elements that are going to be encoded.
     *  @param vs     values of elements that are going to be encoded.
     *  @param newrow vector to append the encoded elements.
     */
    void DoEncodeBlock(std::vector<uint64_t> &xs, std::vector<double> &vs,
                       std::vector<SpmRowElem> &newrow);

    /**
     *  Encode elements from a block when the split_blocks function is active.
     *  @param xs     deltas of elements that are going to be encoded.
     *  @param vs     values of elements that are going to be encoded.
     *  @param newrow vector to append the encoded elements.
     */
    void DoEncodeBlockAlt(std::vector<uint64_t> &xs, std::vector<double> &vs,
                          std::vector<SpmRowElem> &newrow);
    /**
     *  Decode CSX element.
     *  @param elem   element to decode.
     *  @param newrow vector to append the decoded elements.
     */
    void DoDecode(const SpmRowElem *elem, std::vector<SpmRowElem> &newrow);

    /**
     *  Encode the elements of a row.
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
     *  @param rstart first element of row.
     *  @param rend   last element of row.
     *  @param newrow vector to append the decoded elements.
     */
    void DecodeRow(const SpmRowElem *rstart, const SpmRowElem *rend,
	           std::vector<SpmRowElem> &newrow);

    /**
     *  Cut the blocks that pass the max elements limit.
     *  @param stats       stats of the matrix.
     *  @param block_align block's main dimension.
     */ 
    void CutMaxLimit(DeltaRLE::Stats *stats, uint64_t block_align);
    
    /**
     *  Update stats taking into account the split block function.
     *  @param stats       stats of the matrix.
     *  @param size        number of non-zero elementw of initial matrix.
     *  @param block_align block's main dimension.
     */
    void HandleStats(DeltaRLE::Stats *stats, uint64_t size,
                      uint64_t block_align);

    /**
     *  Update the stats of a sub-matrix.
     *  @param spm   object of spm class that describes the sub-matrix elements.
     *  @param xs    deltas of elements examined.
     *  @param stats stats of the sub-matrix so far.
     */
    void UpdateStats(SPM *spm, std::vector<uint64_t> &xs,
                     DeltaRLE::Stats &stats);
    void UpdateStats(std::vector<uint64_t> &xs, DeltaRLE::Stats &stats);
    
    /**
     *  Add stats of a type.
     *  @param stats stats added.
     *  @param type  type which stats are being updated.
     */
    void UpdateStats(SpmIterOrder type, DeltaRLE::Stats stats);
    
    /**
     *  Updates the stats for a block pattern.
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

    //TODO : Ask bkk for comments.
    void ComputeSortSplits();
    void CheckAndSetSorting();
    void DoComputeSortSplitsByRows();
    void DoComputeSortSplitsByNNZ();
    void DoCheckSortByRows();
    void DoCheckSortByNNZ();

    void CheckPropability(double probability) {
        if (probability < 0.0 || probability > 1.0)
            assert(false && "invalid sampling probability");
    }
};

/**
 *  Print stats of a sub-matrix.
 *  @param stats stats of all the available types.
 *  @param spm   object of spm class that describes a sub-matrix.
 *  @param os    output channel.
 */
void DRLE_OutStats(DeltaRLE::Stats &stats, SPM &spm, std::ostream &os);

/**
 *  Keeps data of an encoding sequence.
 */
class Node {
public:
    uint32_t depth;
    std::map <SpmIterOrder, std::set<uint64_t> > deltas_path;
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
     *  @param type type which is ignored.
     */
    void Ignore(SpmIterOrder type);
    
    /**
     *  Copies a node to a new one and inserts an extra type in the end of the
     *  encoding sequence.
     *  @param type   type which is inserted in the end.
     *  @param deltas deltas corresponding to type inserted.
     */
    Node MakeChild(SpmIterOrder type, std::set<uint64_t> deltas);
};

}

#endif

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
