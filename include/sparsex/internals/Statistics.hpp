/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Statistics.hpp
 * \brief SparseX statistics collection framework
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_STATISTICS_HPP
#define SPARSEX_INTERNALS_STATISTICS_HPP

#include <sparsex/internals/Encodings.hpp>
#include <boost/format.hpp>
#include <cassert>
#include <iterator>
#include <map>
#include <ostream>
#include <set>
#include <vector>

using namespace std;

namespace sparsex {
namespace csx {

/**
 *  @brief Statistics data collected by CSX.
 *
 *  The statistics data collected during the preprocessing phase of CSX in order
 *  to determine the best matrix encodings. Three quantities are collected:
 *  
 *  - Number of non-zero elements encoded
 *  - Number of patterns encoded
 *  - Number of delta units generated from a specific encoding.
 *
 *  These values will be used for assigning a score to every possible matrix
 *  encoding and eventually decide on the best option.
 *
 *  @remarks This class is immutable.
 *
 *  @see EncodingManager
 */
class StatsData
{
public:
    /**
     *  Construct a new statistics data instance.
     *
     *  @param nr_encoded number of encoded non-zero elements.
     *  @param nr_patterns number of encoded patterns.
     *  @param nr_deltas number of delta units generated.
     *
     *  \throws This constructor offers the no-throw guarantee.
     */
    StatsData(size_t nr_encoded, size_t nr_patterns, size_t nr_deltas = 0)
        : nr_encoded_(nr_encoded),
          nr_patterns_(nr_patterns),
          nr_deltas_(nr_deltas)
    { }

    /**
     *  Default construct a new statistics data instance.
     *
     *  This is equivalent to StatsData(0, 0, 0).
     */
    StatsData()
        : nr_encoded_(0),
          nr_patterns_(0),
          nr_deltas_(0)
    { }

    /**
     *  Get the number of encoded non-zero elements stored in this instance.
     *
     *  @return the number of encoded non-zero elements stored in this instance.
     */
    size_t GetEncoded() const
    {
        return nr_encoded_;
    }

    /**
     *  Get the number of encoded patterns stored in this instance.
     *
     *  @return the number of encoded patterns stored in this instance.
     */
    size_t GetPatterns() const
    {
        return nr_patterns_;
    }

    /**
     *  Get the number of generated deltas units stored in this instance.
     *
     *  @return the number of generated deltas units stored in this instance.
     */
    size_t GetDeltas() const
    {
        return nr_deltas_;
    }

    /**
     *  Accumulate the data of rhs to this instance.
     *
     *  @param rhs the statistics data to accumulate from.
     */ 
    StatsData &operator+=(const StatsData &rhs)
    {
        nr_encoded_ += rhs.nr_encoded_;
        nr_patterns_ += rhs.nr_patterns_;
        nr_deltas_ += rhs.nr_deltas_;
        return *this;
    }

    /**
     *  Scale this statistics data by a constant factor.
     *
     *  Results are truncated downwards to the nearest integer.
     *
     *  @param rhs the constant factor to multiply with.
     *  @throws std::invalid_argument if rhs is negative.
     */
    StatsData &operator*=(double rhs)
    {
        if (rhs < 0)
            throw invalid_argument("cannot scale with a negative");

        nr_encoded_ *= rhs;
        nr_patterns_ *= rhs;
        nr_deltas_ *= rhs;
        return *this;
    }

    /**
     *  Equality operator.
     *
     *  @param rhs the StatsData to compare with.
     */
    bool operator==(const StatsData &rhs)
    {
        return nr_encoded_ == rhs.nr_encoded_ &&
            nr_patterns_ == rhs.nr_patterns_  &&
            nr_deltas_ == rhs.nr_deltas_;
    }

    /**
     *  Inequality operator.
     *
     *  @param rhs the StatsData to compare with.
     */
    bool operator!=(const StatsData &rhs)
    {
        return !(*this == rhs);
    }

    friend ostream &operator<<(ostream &, const StatsData &);

private:
    size_t nr_encoded_;
    size_t nr_patterns_;
    size_t nr_deltas_;
};

inline ostream &operator<<(ostream &out, const StatsData &data)
{
    out << boost::format("[nz:%u, p:%u, d:%u]") %
        data.nr_encoded_ % data.nr_patterns_ % data.nr_deltas_;
    return out;
}

template<typename T>
class StatsManipulator;

/**
 *  @brief Data structure for maintaining CSX statistics.
 *
 *  This class represents the data structure used for collecting statistics data
 *  of type T. StatsCollection is essentially an accumulation tree comprising
 *  three levels:
 *  
 *  - Instantiation level (leaves): the nodes of this level contain the
 *    statistics data of individual instantiations.
 *   
 *  - Type level: these are intermediate nodes representing a specific type;
 *    children of such nodes are specific instantiations of their type.
 *
 *  - Root node: this level contains an accumulation of all the statistics
 *    inside the collection.
 *
 *  Every level of the collection is associated with a statistics data structure
 *  of type T, which must implement operator+(), operator==() and a default
 *  constructor.
 *
 *  Every time new statistics data is added to the collection, the intermediate
 *  nodes and their associated data are created or updated accordingly.
 *
 *  Although iteration of the stored statistics data is supported, direct access
 *  to the internals of StatsCollection is not permitted whatsoever. However,
 *  there is a need to allow the manipulation of the collected statistics, e.g.,
 *  cutting off instantiations with very low coverage of non-zero elements. This
 *  possible through the use of the so-called "statistics manipulators", which
 *  implement a variation of the visitor design pattern: StatsCollection
 *  iterates over its internal nodes and passes references of their data to the
 *  manipulator. The manipulator is then free to alter this data.
 *
 *  Statistics manipulation is possible at both the instantiation and the type
 *  level. At the type level a manipulator has access to all the stored
 *  instantiations of this type, but not to its associated statistics data
 *  structure; this is updated automatically by the StatsCollection in order to
 *  remain always consistent.
 *
 *  @see StatsCollection::ManipulateStats, StatsManipulator
 */
template<typename T>
class StatsCollection
{
public:
    /**
     *  Instantiation statistics.
     *
     *  An std::map mapping specific instantiations to their associated
     *  statistics data.
     */
    typedef std::map<size_t, T> InstStats;
    
    /**
     *  Construct an empty collection.
     */
    StatsCollection()
        : root_(TypeStatsNode())
    { }

    /**
     *  Get the global (root) statistics of this collection.
     *
     *  These statistics are simply an accumulation of all the individual
     *  statistics stored in the collection.
     *
     *  @return the global statistics data of this collection.
     */
    T GetGlobalStats() const
    {
        return root_.GetData();
    }

    /**
     *  Get the accumulated statistics of an individual encoding type.
     *
     *  @param type the encoding type to check for statistics.
     *  @return the accumulated statistics of the requested encoding type or a
     *          default constructed statistics data type T.
     */
    T GetTypeStats(Encoding::Type type) const
    {
        const TypeStats &type_stats = root_.GetStats();
        typename TypeStats::const_iterator ti = type_stats.find(type);
        if (ti == type_stats.end())
            return T();

        return ti->second.GetData();
    }

    /**
     *  Get the statistics of a specific encoding instantiation.
     *
     *  @param inst the encoding instantiation to check for statistics.
     *  @return the statistics of the requested encoding instantiation or a
     *  default constructed statistics data type T.
     */
    T GetInstStats(const Encoding::Instantiation &inst) const
    {
        const TypeStats &type_stats = root_.GetStats();
        typename TypeStats::const_iterator ti = type_stats.find(inst.first);
        if (ti == type_stats.end())
            return T();

        const InstStats &inst_stats = ti->second.GetStats();
        typename InstStats::const_iterator ii =
            inst_stats.find(inst.second);
        if (ii == inst_stats.end())
            return T();

        return ii->second;
    }

    /**
     *  Append new statistics data or update existing one.
     *
     *  If statistics data already exists for the inst instantiation, then data
     *  will be accumulated to the statistics data of the already stored
     *  instantiation. Otherwise, the new instantiation statistics data will
     *  copied and inserted into the collection.
     *
     *  @param inst the encoding instantiation to append
     *  @param data the statistics data of the appending instantiation.
     */
    void AppendStats(const Encoding::Instantiation &inst, const T &data)
    {
        InstStatsNode new_inst(inst, data);
        root_.AppendNode(inst, new_inst);
    }

    /**
     *  Append or replace statistics data from a different collection.
     *
     *  If the requested type is not found in the other collection, then nothing
     *  happens. Otherwise, this collection's type statistics are copied from
     *  other replacing any data stored under the same type.
     *
     *  This method copies the entire statistics data subtree under type from
     *  the other collection.
     *
     *  @param type the encoding type to copy from the other collection.
     *  @param other the collection to copy from
     */
    void AppendStats(Encoding::Type type, const StatsCollection<T> &other);

    /**
     *  Manipulate statistics using a specific manipulator.
     *
     *  This method iterates all the nodes of this statistics collection in a
     *  post order fashion and calls the appropriate
     *  StatsManipulator::Manipulate method of the passed manipulator depending
     *  on the level accessed:
     *  
     * - At the instantiation level, it is passed a reference to the stored
     *   instantiation data.
     * - At the type level, it is passed a reference to the entire type subtree.
     *   This allows certain manipulators (see BlockSplitter) to have a complete
     *   view of the encoded instantiations of a certain type, taking global
     *   decisions.
     *
     *  The root level is not processed by the manipulator.
     *
     *  If the StatsManipulator::Manipulate method of the passed statistics
     *  manipulator returns true, then intermediate node statistics are
     *  recalculated, in order to remain consistent.
     *
     *  A cleanup step is also performed by this method before returning,
     *  erasing all stale instantiations and/or types, i.e., instantiations
     *  and/or types that were zeroed out by the manipulator.
     *
     *  @remark An instantiation or type is considered zeroed if it is equal to
     *  the default constructed statistics data type T.
     *
     *  @param m the manipulator to which to delegate the manipulation of this
     *           collection statistics data.
     */
    void ManipulateStats(StatsManipulator<T> &m);

    /**
     *  Clear the entire statistics data stored in this collection.
     */
    void Clear()
    {
        root_ = TypeStatsNode();
    }

    // StatsCollection's iterator
    class const_iterator;

    /**
     *  Retrieve a constant iterator at the beginning of this collection.
     *
     *  Iteration is supported only on type nodes.
     *
     *  @return a constant iterator at the beginning of this collection.
     */ 
    const_iterator cbegin() const
    {
        return const_iterator(this);
    }

    /**
     *  Retrieve a constant iterator at the end of this collection.
     *
     *  Iteration is supported only on type nodes.
     *  
     *  @return a constant iterator at the end of this collection.
     */
    const_iterator cend() const
    {
        return const_iterator(this, typename const_iterator::end_tag());
    }

    template<typename U>
    friend ostream &operator<<(ostream &, const StatsCollection<U> &);

private:
    class InstStatsNode
    {
    public:
        InstStatsNode(const Encoding::Instantiation &inst, const T& data)
            : data_(data)
        {
            stats_[inst.second] = data_;
        }

        InstStatsNode() {}

        void AppendData(size_t inst, const T& data)
        {
            if (stats_.find(inst) == stats_.end()) 
                stats_[inst] = data;
            else
                stats_[inst] += data;

            data_ += data;
        }

        T GetData() const
        {
            return data_;
        }

        const InstStats &GetStats() const
        {
            return stats_;
        }

        InstStats &GetStats()
        {
            return stats_;
        }

        void RecalculateStats()
        {
            typename InstStats::iterator ii = stats_.begin();
            typename InstStats::iterator ie = stats_.end();
            data_ = T();
            for (; ii != ie; ++ii)
                data_ += ii->second;
        }

        void Erase(size_t inst)
        {
            stats_.erase(inst);
        }

    private:
        InstStats stats_;
        T data_;
    };

    typedef std::map<Encoding::Type, InstStatsNode> TypeStats;

    class TypeStatsNode
    {
    public:
        TypeStatsNode() { }

        void AppendNode(const Encoding::Instantiation &inst,
                        const InstStatsNode &node)
        {
            if (stats_.find(inst.first) == stats_.end())
                stats_[inst.first] = node;
            else
                stats_[inst.first].AppendData(inst.second, node.GetData());

            data_ += node.GetData();
        }

        T GetData() const
        {
            return data_;
        }

        const TypeStats &GetStats() const
        {
            return stats_;
        }

        TypeStats &GetStats()
        {
            return stats_;
        }

        void RecalculateStats()
        {
            typename TypeStats::iterator ii = stats_.begin();
            typename TypeStats::iterator ie = stats_.end();
            data_ = T();
            for (; ii != ie; ++ii)
                data_ += ii->second.GetData();
        }

        void Erase(Encoding::Type type)
        {
            stats_.erase(type);
        }

    private:
        TypeStats stats_;
        T data_;
    };

    void DeleteTypes(const vector<Encoding::Type> &types)
    {
        for (size_t i = 0; i < types.size(); ++i)
            root_.Erase(types[i]);
    }

    void DeleteInstantiations(const vector<Encoding::Instantiation> &inst)
    {
        for (size_t i = 0; i < inst.size(); ++i)
            root_.GetStats()[inst[i].first].Erase(inst[i].second);
    }

    TypeStatsNode root_;
};

template<typename T>
ostream &operator<<(ostream &out, const StatsCollection<T> &collection)
{
    const typename StatsCollection<T>::TypeStats &type_stats =
        collection.root_.GetStats();
    typename StatsCollection<T>::TypeStats::const_iterator ti =
        type_stats.begin();
    for (; ti != type_stats.end(); ++ti) {
        const typename StatsCollection<T>::InstStats &inst_stats =
            ti->second.GetStats();
        typename StatsCollection<T>::InstStats::const_iterator ii =
            inst_stats.begin();
        Encoding e(ti->first);
        out << e.GetShortName() << ":" << ti->second.GetData() << ": { ";
        for (; ii != inst_stats.end(); ++ii)
            out << ii->first << ":" << ii->second << " ";
        out << "}\n";
    }

    out << "Total: " << collection.root_.GetData();
    return out;
}

/**
 *  @brief Abstract base class of statistics manipulators.
 *
 *  Statistics manipulators are responsible for manipulating the
 *  StatsCollections' stored data. They must implement two methods:
 *  
 *  - One for manipulating an encoding instantiation's statistics data and
 *  - one for manipulating all the instantiations of a specific type.
 *
 *  The latter method is useful for cases when the statistics manipulation
 *  process needs a broad look of all the encoded instantiations of a specific
 *  type (see BlockSplitter).
 *
 *  It is not necessary that both methods do useful work. Any method can simply
 *  return false to denote that it has not altered any data. In this case, the
 *  StatsCollection will not update the statistics data of its intermediate
 *  nodes, thus saving some computation time.
 *
 *  @see StatsCollection::Manipulate
 */
template<typename T>
class StatsManipulator
{
public:
    /**
     *  Manipulate a specific instantiation's data.
     *
     *  @param inst the encoding instantiation whose data will be manipulated.
     *  @param data the statistics data associated with inst.
     *  @return true is data was manipulated, false otherwise.
     */
    virtual bool Manipulate(const Encoding::Instantiation &inst, T &data) = 0;

    /**
     *  Manipulate the entire instantiation subtree of a specific type.
     *
     *  @param type the encoding type whose instantiation data will be
     *              manipulated
     *  @param stats the statistics data of all the instantiations associated
     *               with type.
     *
     *  @return true if the data was manipulated, false otherwise. 
     */
    virtual bool Manipulate(Encoding::Type type,
                            typename StatsCollection<T>::InstStats &stats) = 0;
};

template<typename T>
void StatsCollection<T>::ManipulateStats(StatsManipulator<T> &m)
{
    vector<Encoding::Type> types_to_erase;
    vector<Encoding::Instantiation> inst_to_erase;
    TypeStats &type_stats = root_.GetStats();
    typename TypeStats::iterator ti = type_stats.begin();
    int recalc = 0;
    for (; ti != type_stats.end(); ++ti) {
        InstStats &inst_stats = ti->second.GetStats();
        typename InstStats::iterator ii = inst_stats.begin();
        recalc = 0;
        for (; ii != inst_stats.end(); ++ii) {
            recalc += m.Manipulate(make_pair(ti->first, ii->first), ii->second);
            if (ii->second == T())
                inst_to_erase.push_back(make_pair(ti->first, ii->first));
        }

        recalc += m.Manipulate(ti->first, inst_stats);
        if (recalc)
            ti->second.RecalculateStats();

        if (ti->second.GetData() == T())
            types_to_erase.push_back(ti->first);
    }

    DeleteInstantiations(inst_to_erase);
    DeleteTypes(types_to_erase);
    if (recalc)
        root_.RecalculateStats();
}


template<typename T>
void StatsCollection<T>::AppendStats(Encoding::Type type,
                                     const StatsCollection<T> &other)
{
    TypeStats &this_stats = root_.GetStats();
    const TypeStats &other_stats = other.root_.GetStats();
    typename TypeStats::const_iterator ti = other_stats.find(type);
    if (ti == other_stats.end())
        // type not found in other collection
        return;

    this_stats[type] = ti->second;
    root_.RecalculateStats();
}

/**
 *  @brief A simple manipulator for scaling instantiation data by a constant
 *         factor.
 *
 *  @see StatsManipulator
 *  @see EncodingManager
 */ 
class StatsDataScaler : public StatsManipulator<StatsData>
{
public:
    /**
     *  Construct a new StatsDataScaler.
     *
     *  @param factor the constant factor to scale the statistics.
     */ 
    StatsDataScaler(double factor)
        : factor_(factor)
    {}

    /**
     *  Scale the statistics data of the passed instantiation.
     *
     *  @param inst the encoding instantiation to manipulate
     *  @param data the encoding instantiation's associated data
     *  @return true
     */
    bool Manipulate(const Encoding::Instantiation &inst, StatsData &data)
    {
        data *= factor_;
        return true;
    }

    /**
     *  No manipulation.
     *
     *  @return false
     */
    bool Manipulate(Encoding::Type,
                    typename StatsCollection<StatsData>::InstStats &)
    {
        return false;
    }

private:
    double factor_;
};

/**
 *  @brief A simple statistics manipulator for cutting off instantiation below
 *         a coverage threshold.
 *
 *  @see StatsManipulator
 */
class CoverageFilter : public StatsManipulator<StatsData>
{
public:
    /**
     *  Construct a new CoverageFilter manipulator.
     *
     *  @param nr_nzeros the number of non-zero elements to check coverage
     *          against. This is usually the non-zero elements of a
     *          specific matrix partition (see SparsePartition).
     *  @param min_coverage the minimum coverage. All instantiation below this
     *          threshold will be filtered out.
     *  @param unfiltered an std::set to be filled with the encoding
     *          instantiations that were not filtered.
     */
    CoverageFilter(size_t nr_nzeros, double min_coverage,
                   set<Encoding::Instantiation> &unfiltered)
        : nr_nzeros_(nr_nzeros),
          min_coverage_(min_coverage),
          unfiltered_(unfiltered)
    { }

    /**
     *  Filter out instantiations below this filter's threshold.
     *
     *  @param inst the current instantiation to check.
     *  @param data the statistics data associated with inst. If inst is
     *      filtered, data is replaced by default constructed StatsData.
     *  @return true if instantiation was filtered out, false otherwise.
     *  @see StatsData
     */
    bool Manipulate(const Encoding::Instantiation &inst, StatsData &data)
    {
        bool ret = false;
        double coverage = data.GetEncoded() / (double) nr_nzeros_;
        if (coverage < min_coverage_) {
            data = StatsData();
            ret = true;
        } else {
            unfiltered_.insert(inst);
        }

        return ret;
    }

    /**
     *  No manipulation.
     *
     *  @return false
     */
    bool Manipulate(Encoding::Type type,
                    typename StatsCollection<StatsData>::InstStats &stats)
    {
        return false;
    }

private:
    size_t nr_nzeros_;
    double min_coverage_;
    set<Encoding::Instantiation> &unfiltered_;
};

/**
 *  @brief A statistics manipulator implementing the CSX split block technique.
 *
 *  The idea of the split block technique is to split larger blocks to smaller
 *  ones and merge them with existing dominant blocks in terms of non-zero
 *  coverage, so as to increase the total non-zero coverage of block patterns.
 *
 *  This manipulator searches for the first largest instantiation that exceeds a
 *  minimum threshold in non-zeros coverage and splits all larger blocks to a
 *  sum of this target instantiation and a remainder. This procedure continues
 *  until no more instantiations are left to split.
 *
 *  This manipulator uses the same technique to split also large blocks to a
 *  maximum CSX unit size. This is currently set to 255 in EncodingManager.
 *
 *  @see StatsManipulator
 *  @see EncodingManager
 */
class BlockSplitter : public StatsManipulator<StatsData>
{
public:
    /**
     *  Create a new BlockSplitter instance.
     *
     *  @param max_patt_size maximum CSX unit size allowed.
     *  @param nr_nzeros number of non-zero elements to the check the coverage
     *      against
     *  @param min_coverage the minimum coverage threshold
     */
    BlockSplitter(size_t max_patt_size, size_t nr_nzeros, double min_coverage)
        : max_patt_size_(max_patt_size),
          nr_nzeros_(nr_nzeros),
          min_coverage_(min_coverage)
    { }

    /**
     *  No manipulation happens at this level.
     *
     *  @return false
     */
    bool Manipulate(const Encoding::Instantiation &inst, StatsData &data)
    {
        return false;
    }

    /**
     *  Implements the CSX split block technique.
     *
     *  @param type the block type to apply the split block technique upon. If
     *      type does not refer to a block type, this method does not touch
     *      the statistics.
     *  @param stats the statistics data of the instantiations associated with
     *      type
     *  @return true if statistics manipulation occurred, false otherwise.
     */
    bool Manipulate(Encoding::Type type,
                    typename StatsCollection<StatsData>::InstStats &stats);

private:
    size_t max_patt_size_;
    size_t nr_nzeros_;
    double min_coverage_;
};

// StatsCollection iterator
template<typename T>
class StatsCollection<T>::const_iterator :
        public iterator<forward_iterator_tag, pair<Encoding::Type, T> >
{
public:
    struct end_tag {};
    const_iterator(const StatsCollection<T> *stats)
        : stats_(stats),
          curr_(stats->root_.GetStats().begin())
    { }

    const_iterator(const StatsCollection<T> *stats, const end_tag &)
        : stats_(stats),
          curr_(stats->root_.GetStats().end())
    { }

    pair<Encoding::Type, T> operator*() const
    {
        return make_pair(curr_->first, curr_->second.GetData());
    }

    void operator++()
    {
        ++curr_;
    }

    bool operator==(const const_iterator &rhs) const
    {
        return stats_ == rhs.stats_ && curr_ == rhs.curr_;
    }

    bool operator!=(const const_iterator &rhs) const
    {
        return !(*this == rhs);
    }

private:
    const StatsCollection<T> *stats_;
    typename StatsCollection<T>::TypeStats::const_iterator curr_;
};

} // end of namespace csx
} // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_STATISTICS_HPP
