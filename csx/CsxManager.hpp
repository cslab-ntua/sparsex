/* -*- C++ -*-
 *
 * CsxManager.hpp -- The CSX Manager interface.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_MANAGER_HPP
#define CSX_MANAGER_HPP

#include "Allocators.hpp"
#include "Csx.hpp"
#include "CtlBuilder.hpp"
#include "CtlUtil.hpp"
#include "Element.hpp"
#include "SparsePartition.hpp"
#include "Delta.hpp"
#include "SpmMt.hpp"
#include "Logger.hpp"

#ifdef SPM_NUMA
#   include <numa.h>
#   include "numa_util.h"
#endif

using namespace std;

namespace csx {

/**
 *  This class is responsible for transforming the matrix from our internal
 *  representation (SPM object) to the final CSX form.
 *
 *  @see SPM
 */
template<typename IndexType, typename ValueType>
class CsxManager
{
public:
    /**
     *  Pattern information that is essential for the construction of CSX.
     */
    struct PatInfo {
        PatInfo(uint8_t flag_, uint64_t npatterns_, uint64_t nr_): 
                flag(flag_), npatterns(npatterns_), nr(nr_) {}
        PatInfo(): flag(0), npatterns(0), nr(0) {}

        uint8_t flag;       // A unique CSX ID assigned to this pattern.
        uint64_t npatterns; // Number of patterns of this kind.
        uint64_t nr;        // Number of non-zero elements of this pattern.
    };

    typedef std::map<long,PatInfo> PatMap;
    PatMap patterns;        // Patterns found in matrix.

    CsxManager(SparsePartition<IndexType, ValueType> *spm)
        : spm_(spm),
          spm_sym_(0),
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
          ctl_builder_(CtlBuilder(6*spm_->GetNrNonzeros())),
          last_col_(0), empty_rows_(0) {}

    CsxManager(SparsePartitionSym<IndexType, ValueType> *spm_sym)
        : spm_(0),
          spm_sym_(spm_sym),
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
          ctl_builder_(
              CtlBuilder(6*spm_sym_->GetLowerMatrix()->GetNrNonzeros())),
          last_col_(0), empty_rows_(0) {}

    /**
     *  Get a unique CSX ID for the pattern with SPM ID <tt>pattern_id</tt> and
     *  updates statistics for this pattern.
     *
     *  @param pattern_id the pattern ID in the SPM format.
     *  @param nnz        number of elements included in the specific pattern.
     *  @return           the CSX pattern ID for the specified SPM pattern.
     *  @see SPM
     */
    uint8_t GetFlag(unsigned long pattern_id, size_t nnz);

    /**
     *  Transform the matrix owned by this manager into CSX form.
     *  
     *  @param symmetric determines if the matrix is symmetric (only the lower
     *                   triangle half must be encoded) or not.
     *  @return          a handle to the newly created CSX matrix or to CSX-Sym
     *                   lower triangle half part of the matrix.
     */
    csx_t<ValueType> *MakeCsx(bool symmetric);
    
    /**
     *  Transform the matrix owned by this manager into CSX form.
     *
     *  @return a handle to the newly created CSX-Sym matrix.
     */
    csx_sym_t<ValueType> *MakeCsxSym();

    /**
     *  Checks whether row jumps exist in the matrix to be encoded in CSX
     *  format.
     *
     *  @return <tt>true</tt> if row jumps exist.
     */
    bool HasRowJmps()
    {
        return row_jmps_;
    }

    bool HasFullColumnIndices()
    {
        return full_column_indices_;
    }

    void SetFullColumnIndices(bool val)
    {
        full_column_indices_ = val;
    }

private:
    void AddMappings(long *map);
    /**
     *  Transform a row of the matrix into CSX form.
     *
     *  @rstart first element of the row.
     *  @rend   last element of the row.
     */
    void DoRow(
        typename SparsePartition<IndexType, ValueType>::iterator &rstart,
        typename SparsePartition<IndexType, ValueType>::iterator &rend,
        IndexType row);
    void DoSymRow(
        typename SparsePartition<IndexType, ValueType>::iterator &rstart,
        typename SparsePartition<IndexType, ValueType>::iterator &rend);

    /**
     *  Set flags that concern change of row.
     *
     *  @return nr bit + ujmp (if ujmp = 0, no empty rows then)
     */
    pair<bool, size_t> UpdateNewRow();

    /**
     *  Set all the flags of pattern info for elements that do not adhere to
     *  any type of pattern.
     *
     *  @param xs vector with the columns of elements.
     */
    void AddXs(std::vector<IndexType> &xs);

    /**
     *  Set all the flags of pattern info for elements that adhere to a type of
     *  pattern.
     *
     *  @param elem elements of current pattern.
     */
    void AddPattern(const Element<IndexType, ValueType> &elem);

    /**
     *  @param xs  elements found before pattern.
     */
    void PreparePat(std::vector<IndexType> &xs,
                    const Element<IndexType, ValueType> &elem);

    /**
     *  Updates the span of a row.
     *
     *  @param elem     An element representing a pattern.
     */
    void UpdateRowSpan(const Element<IndexType, ValueType> &elem);

    SparsePartition<IndexType, ValueType> *spm_;
    SparsePartitionSym<IndexType, ValueType> *spm_sym_;
    uint8_t flag_avail_;     // Current available flags for pattern id mapping.
    bool row_jmps_;          // Whether or not row jumps included.
    bool full_column_indices_;  // use full 32-bit indices instead of deltas.
                                
    ValueType *values_;
    row_info_t *rows_info_; 
    size_t span_;
    size_t values_idx_;

    CtlBuilder ctl_builder_;
    IndexType last_col_;
    bool new_row_;          // Marker of new_row.
    uint64_t empty_rows_;   // Number of empty rows since last non-empty row.
    size_t curr_row_;
};

}   // end of csx namespace

/* Helper functions */
template<typename IterT, typename ValT>
void DeltaEncode(IterT start, IterT end, ValT &x0)
{
    IterT i;
    ValT prev, tmp;

    prev = x0;
    for (i = start; i != end; ++i) {
        tmp = *i;
        *i -= prev;
        prev = tmp;
    }
}

template<typename S, typename D>
void Copy(D *dst, S *src, long nr_items)
{
    for (long i = 0; i < nr_items; i++){
        dst[i] = static_cast<D>(src[i]);
    }
}

namespace csx {

/* 
 * Implementation of class CsxManager
 */
template<typename IndexType, typename ValueType>
uint8_t CsxManager<IndexType, ValueType>::GetFlag(unsigned long pattern_id,
                                                  size_t nnz)
{
    typename CsxManager::PatMap::iterator pi;
    uint8_t ret;

    pi = this->patterns.find(pattern_id);
    if (pi == this->patterns.end()) {
        ret = flag_avail_++;
        assert(ret <= CTL_PATTERNS_MAX && "too many patterns");
        typename CsxManager<IndexType, ValueType>::PatInfo pat_info(ret, 
                                                                    1, nnz);
        this->patterns[pattern_id] = pat_info;
    } else {
        ret = pi->second.flag;
        pi->second.npatterns++;
        pi->second.nr += nnz;
    }

    return ret;
}

template<typename IndexType, typename ValueType>
csx_sym_t<ValueType> *CsxManager<IndexType, ValueType>::MakeCsxSym()
{
    csx_sym_t<ValueType> *csx;
    double *diagonal = spm_sym_->GetDiagonal();
    uint64_t diagonal_size = spm_sym_->GetDiagonalSize();
#ifdef SPM_NUMA
    NumaAllocator &numa_alloc = NumaAllocator::GetInstance();
#endif
    
    spm_ = spm_sym_->GetLowerMatrix();

#ifdef SPM_NUMA
    int cpu = sched_getcpu();
    if (cpu < 0) {
        LOG_ERROR << "sched_getcpu() failed " << strerror(errno);
        exit(1);
    }

    int node = numa_node_of_cpu(cpu);
    if (node < 0) {
        LOG_ERROR << "numa_node_of_cpu() failed " << strerror(errno);
        exit(1);
    }

    csx = new (numa_alloc, node) csx_sym_t<ValueType>;
    csx->dvalues = new (numa_alloc, node) ValueType[diagonal_size];
#else  
    csx = new csx_sym_t<ValueType>;
    csx->dvalues = new ValueType[diagonal_size];
#endif

    for (uint64_t i = 0; i < diagonal_size; i++)
        csx->dvalues[i] = diagonal[i];
    csx->lower_matrix = MakeCsx(true);
    
    return csx;
}

template<typename IndexType, typename ValueType>
csx_t<ValueType> *CsxManager<IndexType, ValueType>::MakeCsx(bool symmetric)
{
    csx_t<ValueType> *csx;

#ifdef SPM_NUMA
    NumaAllocator &numa_alloc = NumaAllocator::GetInstance();
    int cpu = sched_getcpu();
    if (cpu < 0) {
        LOG_ERROR << "sched_getcpu() failed " << strerror(errno);
        exit(1);
    }

    int node = numa_node_of_cpu(cpu);
    if (node < 0) {
        LOG_ERROR << "numa_node_of_cpu() failed " << strerror(errno);
        exit(1);
    }

    csx = new (numa_alloc, node) csx_t<ValueType>;
    values_ = new (numa_alloc, node) ValueType[spm_->GetNrNonzeros()];
    rows_info_ = new (numa_alloc, node) row_info_t[spm_->GetNrNonzeros()];
#else    
    csx = new csx_t<ValueType>;
    values_ = new ValueType[spm_->GetNrNonzeros()];
    rows_info_ = new row_info_t[spm_->GetNrNonzeros()];
#endif  // SPM_NUMA
    if (!csx || !values_ || !rows_info_) {
        LOG_ERROR << "malloc failed\n";
        exit(1);
    }

    // Be greedy with the initial capacity (equal to CSR col_ind size) 
    // to avoid realloc()'s.
    csx->nnz = spm_->GetNrNonzeros();
    csx->nrows = spm_->GetNrRows();
    csx->ncols = spm_->GetNrCols();
    csx->row_start = spm_->GetRowStart();
    values_idx_ = 0;
    new_row_ = false;   // Do not mark first row.

    if (!symmetric) {
        curr_row_ = 0;
        for (size_t i = 0; i < spm_->GetRowptrSize() - 1; ++i, ++curr_row_) {
            typename SparsePartition<IndexType, ValueType>::iterator rbegin =
                spm_->begin(i);
            typename SparsePartition<IndexType, ValueType>::iterator rend =
                spm_->end(i);

            LOG_DEBUG << "MakeCsx(): row: " << i << "\n";
    
            if (rbegin == rend) {		// Check if row is empty.
                LOG_DEBUG << "MakeCsx(): row is empty\n";
                if (new_row_ == false) {
                    rows_info_[i].rowptr = 0;
                    new_row_ = true;	// In case the first row is empty.
                } else {
                    empty_rows_++;
                    rows_info_[i].rowptr = rows_info_[i-1].rowptr;
                }
                rows_info_[i].valptr = 0;
                rows_info_[i].span = 0;
                continue;
            }

            if (i > 0)
                rows_info_[i].rowptr = ctl_builder_.GetCtlSize();
            else 
                rows_info_[i].rowptr = 0;

            rows_info_[i].valptr = values_idx_;
            DoRow(rbegin, rend, i);
            rows_info_[i].span = span_;
            new_row_ = true;
        }

        for (size_t i = spm_->GetRowptrSize() - 1;
             i < (size_t) spm_->GetNrRows(); i++) {
                rows_info_[i].valptr = 0;
                rows_info_[i].rowptr = rows_info_[i-1].rowptr;
                rows_info_[i].span = 0;
        }
    } else {
        for (size_t i = 0; i < spm_->GetRowptrSize() - 1; i++) {
            typename SparsePartition<IndexType, ValueType>::iterator rbegin =
                spm_->begin(i);
            typename SparsePartition<IndexType, ValueType>::iterator rend =
                spm_->end(i);

            LOG_DEBUG << "MakeCsx(): row: " << i << "\n";
    
            if (rbegin == rend){		// Check if row is empty.
                LOG_DEBUG << "MakeCsx(): row is empty\n";
    
                if (new_row_ == false) {
                    rows_info_[i].rowptr = 0;
                    new_row_ = true;	// In case the first row is empty.
                } else {
                    rows_info_[i].rowptr = rows_info_[i-1].rowptr;
                    empty_rows_++;
                }
                rows_info_[i].valptr = 0;
                rows_info_[i].span = 0;
                continue;
            }
    
            if (i > 0)
                rows_info_[i].rowptr = ctl_builder_.GetCtlSize();
            else 
                rows_info_[i].rowptr = 0;
            rows_info_[i].valptr = values_idx_;
            DoSymRow(rbegin, rend);
            rows_info_[i].span = span_;
            new_row_ = true;
        }

        for (size_t i = spm_->GetRowptrSize() - 1;
             i < (size_t) spm_->GetNrRows(); i++) {
                rows_info_[i].valptr = 0;
                rows_info_[i].rowptr = rows_info_[i-1].rowptr;
                rows_info_[i].span = 0;
        }
    }

#ifdef DEBUG_MODE
    LOG_DEBUG << "values_\n";
    for (size_t i = 0; i < spm_->GetNrNonzeros(); ++i)
        LOG_DEBUG << values_[i] << "\n";
#endif

    csx->row_jumps = row_jmps_;
    csx->ctl_size = ctl_builder_.GetCtlSize();
    csx->ctl = ctl_builder_.Finalize();
    assert(values_idx_ == spm_->GetNrNonzeros());
    csx->values = values_;
    values_ = NULL;
    values_idx_ = 0;
    csx->rows_info = rows_info_;
    rows_info_ = NULL;
    AddMappings(csx->id_map);
    return csx;
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::AddMappings(long *map)
{
    typename CsxManager::PatMap::iterator pi = patterns.begin();
    typename CsxManager::PatMap::iterator pi_end = patterns.end();

    for (;pi != pi_end; pi++) {
        map[pi->second.flag] = pi->first;
    }
    map[patterns.size()] = -1;
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::UpdateRowSpan(
    const Element<IndexType, ValueType> &elem)
{
    assert(elem.IsPattern());
    const Encoding::Instantiation &inst = elem.GetInstantiation();
    size_t size = elem.GetSize();
    Encoding::Type type = inst.first;
    size_t delta = inst.second;

    Encoding e(type);
    size_t span;
    switch (type) {
    case Encoding::Vertical:
    case Encoding::Diagonal:
    case Encoding::AntiDiagonal:
        span = (size - 1) * delta;
        break;
    case Encoding::BlockRowMin ... Encoding::BlockRowMax:
        span = type - Encoding::BlockRowMin;
        break;
    case Encoding::BlockColMin ... Encoding::BlockColMax:
        span = size / e.GetBlockAlignment() - 1;
        break;
    default:
        span = 0;
    }

    if (span > span_)
        span_ = span;
}

/*
 *  Ctl Rules
 *  1. Each unit leaves the x index at the last element it calculated on the
 *     current row.
 *  2. Size is the number of elements that will be calculated.
 */
template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
DoRow(typename SparsePartition<IndexType, ValueType>::iterator &rbegin,
      typename SparsePartition<IndexType, ValueType>::iterator &rend,
      IndexType row)
{
    std::vector<IndexType> xs;

    span_ = 0;
    last_col_ = 1;
    typename SparsePartition<IndexType, ValueType>::iterator &ri = rbegin;
    for (; ri != rend; ++ri) {
        LOG_DEBUG << "\t" << *ri << "\n";

        // Check if this element contains a pattern.
        if ((*ri).IsPattern()) {
            UpdateRowSpan(*ri);
            PreparePat(xs, *ri);
            assert(xs.size() == 0);
            AddPattern(*ri);
            const ValueType *elem_vals = &(*ri).GetValues();
            std::copy(elem_vals, elem_vals + (*ri).GetSize(),
                      values_ + values_idx_);
            values_idx_ += (*ri).GetSize();
            continue;
        }

        // Check if we exceeded the maximum size for a unit.
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX)
             AddXs(xs);

        xs.push_back((*ri).GetCol());
        values_[values_idx_++] = (*ri).GetValue();
    }

    if (xs.size() > 0)
        AddXs(xs);
}

/*
 *  Ctl Rules
 *  1. Each unit leaves the x index at the last element it calculated on the
 *     current row.
 *  2. Size is the number of elements that will be calculated.
 */
template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
DoSymRow(typename SparsePartition<IndexType, ValueType>::iterator &rstart,
         typename SparsePartition<IndexType, ValueType>::iterator &rend)
{
    std::vector<IndexType> xs;
    typename SparsePartition<IndexType, ValueType>::iterator &ri = rstart;
    span_ = 0;

    last_col_ = 1;
    for ( ; ri != rend && (*ri).GetCol() < spm_->GetRowStart() + 1; ++ri) {
        LOG_DEBUG << "\t" << *ri << "\n";

        // Check if this element contains a pattern.
        if ((*ri).IsPattern()) {
            UpdateRowSpan(*ri);
            PreparePat(xs, *ri);
            assert(xs.size() == 0);
            AddPattern(*ri);
            const ValueType *elem_vals = &(*ri).GetValues();
            std::copy(elem_vals, elem_vals + (*ri).GetSize(),
                      values_ + values_idx_);
            values_idx_ += (*ri).GetSize();
            continue;
        }

        // Check if we exceeded the maximum size for a unit.
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX)
             AddXs(xs);

        xs.push_back((*ri).GetCol());
        values_[values_idx_++] = (*ri).GetValue();
    }

    if (xs.size() > 0)
        AddXs(xs);

    // FIXME: this is completely the same as before!
    for ( ; ri != rend; ++ri) {
        LOG_DEBUG << "\t" << *ri << "\n";

        // Check if this element contains a pattern.
        if ((*ri).IsPattern()) {
            UpdateRowSpan(*ri);
            PreparePat(xs, *ri);
            assert(xs.size() == 0);
            AddPattern(*ri);
            const ValueType *elem_vals = &(*ri).GetValues();
            std::copy(elem_vals, elem_vals + (*ri).GetSize(),
                      values_ + values_idx_);
            values_idx_ += (*ri).GetSize();
            continue;
        }

        // Check if we exceeded the maximum size for a unit.
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX)
             AddXs(xs);

        xs.push_back((*ri).GetCol());
        values_[values_idx_++] = (*ri).GetValue();
    }

    if (xs.size() > 0)
        AddXs(xs);
}

template<typename IndexType, typename ValueType>
pair<bool, size_t>
CsxManager<IndexType, ValueType>::UpdateNewRow()
{
    bool nr = false;
    size_t rowjmp = 0;

    if (new_row_) {
        nr = true;
        new_row_ = false;
        if (empty_rows_ != 0) {
            rowjmp = empty_rows_ + 1;
            empty_rows_ = 0;
            row_jmps_ = true;
        }
    }

    return make_pair(nr, rowjmp);
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::AddXs(std::vector<IndexType> &xs)
{
    typename std::vector<IndexType>::iterator vi;

    size_t xs_size = xs.size();
    IndexType last_col = xs[xs_size-1];
    IndexType x_start = xs[0];

    // Do delta encoding
    DeltaEncode(xs.begin(), xs.end(), last_col_);
    last_col_ = last_col;

    // Calculate the delta's size and the pattern id
    IndexType max = 0;
    if (xs_size > 1) {
        vi = xs.begin();
        std::advance(vi, 1);
        max = *(std::max_element(vi, xs.end()));
    }

    size_t delta_bytes = GetDeltaSize(max);
    unsigned long patt_id = GetPatternId(delta_bytes);
    
    // Variables ctls_size, ctl_flags are not valid after this call
    pair<bool, size_t> newrow_info = UpdateNewRow();
    IndexType ucol;
    if (full_column_indices_)
        ucol = x_start - 1;
    else
        ucol = xs[0];

    LOG_DEBUG << "AddXs() ucol " << ucol << "\n";
    LOG_DEBUG << "last_col_: " << last_col_ << "\n";

    assert(xs_size > 0 && xs_size <= CTL_SIZE_MAX);
    ctl_builder_.AppendCtlHead(newrow_info.first, newrow_info.second,
                               GetFlag(patt_id, xs_size),
                               xs_size, ucol, sizeof(IndexType),
                               full_column_indices_);

    // Add deltas (if needed)
    for (size_t i = 1; i < xs_size; ++i)
        ctl_builder_.AppendFixedInt(xs[i], delta_bytes);

    xs.clear();
    return;
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
AddPattern(const Element<IndexType, ValueType> &elem)
{
    unsigned long patt_id = GetPatternId(elem);
    size_t patt_size = elem.GetSize();
    pair<bool, size_t> newrow_info = UpdateNewRow();
    IndexType ucol;
    if (full_column_indices_)
        ucol = elem.GetCol() - 1;
    else
        ucol = elem.GetCol() - last_col_;

    assert(patt_size <= CTL_SIZE_MAX);
    ctl_builder_.AppendCtlHead(newrow_info.first, newrow_info.second,
                               GetFlag(patt_id, patt_size),
                               patt_size, ucol, sizeof(IndexType),
                               full_column_indices_);

    LOG_DEBUG << "AddPattern() ucol " << ucol << "\n";
    last_col_ = GetLastCol(elem, spm_->GetType());
    LOG_DEBUG << "last_col_: " << last_col_ << "\n";
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
PreparePat(std::vector<IndexType> &xs,
           const Element<IndexType, ValueType> &elem)
{
    if (xs.size() != 0)
        AddXs(xs);
}

} // end of csx namespace

#endif  // CSX_MANAGER_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
