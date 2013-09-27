/* -*- C++ -*-
 *
 * csx.h -- The CSX Manager interface
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_H__
#define CSX_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "ctl_ll.h"
#include <stdio.h>

#ifdef __cplusplus
}
#endif

///< CSX matrix format
typedef struct {
    size_t rowptr;    /* rowptr is the index in csx->ctl of
                         the first element of row i */
    size_t valptr;    /* valptr is the index in csx->values of
                         the first element of row i */
    size_t span;
} row_info_t;

//padding issue with size_t???
typedef struct {
    size_t nnz, ncols, nrows, ctl_size, row_start;
    uint8_t row_jumps;
    double *values;
    uint8_t *ctl;
    long id_map[CTL_PATTERNS_MAX];
    row_info_t *rows_info; 
} csx_double_t;

typedef struct {
    csx_double_t *lower_matrix;
    double *dvalues;
} csx_double_sym_t;

#ifdef __cplusplus

template<typename ValueType>
struct csx_t {
    size_t nnz, ncols, nrows, ctl_size, row_start;
    uint8_t row_jumps;
    ValueType *values;
    uint8_t *ctl;
    long id_map[CTL_PATTERNS_MAX];
    row_info_t *rows_info;
};

template<typename ValueType>
struct csx_sym_t {
    csx_t<ValueType> *lower_matrix;
    ValueType *dvalues;
};

#endif

#include "dynarray.h"

#define PID_DELTA_BASE 0

#ifdef __cplusplus  // This file is also included from csx_llvm_tmpl.c.

#include "sparse_partition.h"
#include "delta.h"
#include "spm_mt.h"
#include "logger.hpp"

#ifdef SPM_NUMA
#   include <numa.h>
#   include "numa_util.h"
#   define DYNARRAY_CREATE  dynarray_create_numa
#else
#   define DYNARRAY_CREATE  dynarray_create
#endif

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
          ctl_da_(NULL),
          last_col_(0), empty_rows_(0) {}

    CsxManager(SparsePartitionSym<IndexType, ValueType> *spm_sym)
        : spm_(0),
          spm_sym_(spm_sym),
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
          ctl_da_(NULL),
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
    uint8_t GetFlag(long pattern_id, IndexType nnz);

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
    void DoRow(const Elem<IndexType, ValueType> *rstart,
               const Elem<IndexType, ValueType> *rend, IndexType row);
    void DoSymRow(const Elem<IndexType, ValueType> *rstart,
                  const Elem<IndexType, ValueType> *rend);

    /**
     *  Set flags that concern change of row.
     *
     *  @flags flags that are going to be updated.
     */
    void UpdateNewRow(uint8_t *flags);

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
    void AddPattern(const Elem<IndexType, ValueType> &elem);

    /**
     *  @param xs  elements found before pattern.
     */
    void PreparePat(std::vector<IndexType> &xs,
                    const Elem<IndexType, ValueType> &elem);

    /**
     *  Updates the span of a row.
     *
     *  @param pat  the pattern of an element in the row.
     */
    void UpdateRowSpan(const DeltaRLE *pat);

    SparsePartition<IndexType, ValueType> *spm_;
    SparsePartitionSym<IndexType, ValueType> *spm_sym_;
    uint8_t flag_avail_;     // Current available flags for pattern id mapping.
    bool row_jmps_;          // Whether or not row jumps included.
    bool full_column_indices_;  // use full 32-bit indices instead of deltas.
                                
    ValueType *values_;
    row_info_t *rows_info_; 
    size_t span_;
    IndexType values_idx_;

    dynarray_t *ctl_da_;
    IndexType last_col_;
    bool new_row_;          // Marker of new_row.
    uint64_t empty_rows_;   // Number of empty rows since last non-empty row.
};

}   // end of csx namespace

/* Helper functions */
template<typename IterT, typename ValT>
void DeltaEncode(IterT start, IterT end, ValT &x0)
{
    IterT i;
    ValT prev, tmp;

    prev = x0;
    for (i = start; i != end; ++i){
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
uint8_t CsxManager<IndexType, ValueType>::GetFlag(long pattern_id,
                                                  IndexType nnz)
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

    csx = (csx_sym_t<ValueType> *) numa_alloc_onnode
        (sizeof(csx_sym_t<ValueType>), node);
    csx->dvalues = (ValueType *) numa_alloc_onnode(diagonal_size *
                                                   sizeof(ValueType), node);
#else  
    csx = (csx_sym_t<ValueType> *) malloc(sizeof(csx_sym_t<ValueType>));
    csx->dvalues = (ValueType *) malloc(diagonal_size * sizeof(ValueType));
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

    csx = (csx_t<ValueType> *) alloc_onnode
        (sizeof(csx_t<ValueType>), node);
    values_ = (ValueType *) alloc_onnode
        (sizeof(ValueType)*spm_->GetNrNonzeros(), node);
    rows_info_ = (row_info_t *) alloc_onnode
        (sizeof(row_info_t)*spm_->GetNrRows(), node);
#else    
    csx = (csx_t<ValueType> *) malloc(sizeof(csx_t<ValueType>));
    values_ = (ValueType *) malloc(sizeof(ValueType)*spm_->GetNrNonzeros());
    rows_info_ = (row_info_t *) malloc(sizeof(row_info_t)*spm_->GetNrRows());
#endif  // SPM_NUMA
    if (!csx || !values_ || !rows_info_) {
        LOG_ERROR << "malloc failed\n";
        exit(1);
    }

    // Be greedy with the initial capacity (equal to CSR col_ind size) 
    // to avoid realloc()'s.
    ctl_da_ = DYNARRAY_CREATE(sizeof(uint8_t), 512, 6*spm_->GetNrNonzeros());
    csx->nnz = spm_->GetNrNonzeros();
    csx->nrows = spm_->GetNrRows();
    csx->ncols = spm_->GetNrCols();
    csx->row_start = spm_->GetRowStart();
    values_idx_ = 0;
    new_row_ = false;		        // Do not mark first row.

    if (!symmetric) {
        for (size_t i = 0; i < spm_->GetRowptrSize() - 1; i++) {
            const Elem<IndexType, ValueType> *rbegin, *rend;

            rbegin = spm_->RowBegin(i);
            rend = spm_->RowEnd(i);

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
                rows_info_[i].rowptr = (IndexType) dynarray_get_nextidx(ctl_da_);
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
    } 
    else {
        for (size_t i = 0; i < spm_->GetRowptrSize() - 1; i++) {
            const Elem<IndexType, ValueType> *rbegin, *rend;
    
            rbegin = spm_->RowBegin(i);
            rend = spm_->RowEnd(i);
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
                rows_info_[i].rowptr = (IndexType) dynarray_get_nextidx(ctl_da_);
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

    csx->row_jumps = row_jmps_;
    csx->ctl_size = dynarray_size(ctl_da_);
    csx->ctl = (uint8_t *) dynarray_destroy(ctl_da_);
    ctl_da_ = NULL;
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
void CsxManager<IndexType, ValueType>::UpdateRowSpan(const DeltaRLE *pat)
{
    Encoding::Type type = pat->GetType();
    long size = pat->GetSize();
    uint32_t delta = pat->GetDelta();
    size_t span;

    switch (type) {
    case Encoding::Vertical:
    case Encoding::Diagonal:
    case Encoding::AntiDiagonal:
        span = (size - 1) * delta;
        break;
    case Encoding::BlockRowMin ... Encoding::BlockRowMax:
        span = type - BLOCK_TYPE_START - 1;
        break;
    case Encoding::BlockColMin ... Encoding::BlockColMax:
        span = pat->GetOtherDim() - 1;
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
DoRow(const Elem<IndexType, ValueType> *rbegin,
      const Elem<IndexType, ValueType> *rend, IndexType row)
{
    std::vector<IndexType> xs;
    span_ = 0;

    last_col_ = 1;
    for (const Elem<IndexType, ValueType> *spm_elem = rbegin;
         spm_elem < rend; spm_elem++) {
        LOG_DEBUG << "\t" << *spm_elem << "\n";

        // Check if this element contains a pattern.
        if (spm_elem->pattern != NULL) {
            UpdateRowSpan(spm_elem->pattern);
            PreparePat(xs, *spm_elem);
            assert(xs.size() == 0);
            AddPattern(*spm_elem);
            for (long i = 0; i < spm_elem->pattern->GetSize(); i++)
                values_[values_idx_++] = spm_elem->vals[i];

            continue;
        }

        // Check if we exceeded the maximum size for a unit.
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX)
             AddXs(xs);

        xs.push_back(spm_elem->col);
        values_[values_idx_++] = spm_elem->val;
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
DoSymRow(const Elem<IndexType, ValueType> *rbegin,
         const Elem<IndexType, ValueType> *rend)
{
    std::vector<IndexType> xs;
    const Elem<IndexType, ValueType> *spm_elem = rbegin;
    span_ = 0;

    last_col_ = 1;
    for ( ; spm_elem < rend && spm_elem->col < spm_->GetRowStart() + 1; 
         spm_elem++) {
        LOG_DEBUG << "\t" << *spm_elem << "\n";

        // Check if this element contains a pattern.
        if (spm_elem->pattern != NULL) {
            UpdateRowSpan(spm_elem->pattern);
            PreparePat(xs, *spm_elem);
            assert(xs.size() == 0);
            AddPattern(*spm_elem);
            for (long i=0; i < spm_elem->pattern->GetSize(); i++)
                values_[values_idx_++] = spm_elem->vals[i];

            continue;
        }

        // Check if we exceeded the maximum size for a unit.
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX)
             AddXs(xs);

        xs.push_back(spm_elem->col);
        values_[values_idx_++] = spm_elem->val;
    }

    if (xs.size() > 0)
        AddXs(xs);

    for ( ; spm_elem < rend; spm_elem++) {
        LOG_DEBUG << "\t" << *spm_elem << "\n";

        // Check if this element contains a pattern.
        if (spm_elem->pattern != NULL) {
            UpdateRowSpan(spm_elem->pattern);
            PreparePat(xs, *spm_elem);
            assert(xs.size() == 0);
            AddPattern(*spm_elem);
            for (long i=0; i < spm_elem->pattern->GetSize(); i++)
                values_[values_idx_++] = spm_elem->vals[i];

            continue;
        }

        // Check if we exceeded the maximum size for a unit.
        assert(xs.size() <= CTL_SIZE_MAX);
        if (xs.size() == CTL_SIZE_MAX)
             AddXs(xs);

        xs.push_back(spm_elem->col);
        values_[values_idx_++] = spm_elem->val;
    }

    if (xs.size() > 0)
        AddXs(xs);
}

// Note that this function may allocate space in ctl_da.
template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::UpdateNewRow(uint8_t *flags)
{
	if (!new_row_)
		return;

	set_bit(flags, CTL_NR_BIT);
	new_row_ = false;
	if (empty_rows_ != 0){
		set_bit(flags, CTL_RJMP_BIT);
		da_put_ul(ctl_da_, empty_rows_ + 1);
		empty_rows_ = 0;
		row_jmps_ = true;
	}
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::AddXs(std::vector<IndexType> &xs)
{
    uint8_t *ctl_flags, *ctl_size;
    long pat_id, xs_size, delta_bytes;
    IndexType last_col, max;
    DeltaSize delta_size;
    typename std::vector<IndexType>::iterator vi;
    void *dst;

    // Do delta encoding.
    xs_size = xs.size();
    last_col = xs[xs_size-1];
    IndexType x_start = xs[0];
    DeltaEncode(xs.begin(), xs.end(), last_col_);
    last_col_ = last_col;

    // Calculate the delta's size and the pattern id.
    max = 0;
    if (xs_size > 1) {
        vi = xs.begin();
        std::advance(vi, 1);
        max = *(std::max_element(vi, xs.end()));
    }
    delta_size =  getDeltaSize(max);
    pat_id = (8<<delta_size) + PID_DELTA_BASE;

    // Set flags.
    ctl_flags = (uint8_t *) dynarray_alloc_nr(ctl_da_, 2);
    *ctl_flags = GetFlag(PID_DELTA_BASE + pat_id, xs_size);

    // Set size.
    ctl_size = ctl_flags + 1;
    assert( (xs_size > 0) && (xs_size <= CTL_SIZE_MAX));
    *ctl_size = xs_size;

    // Variables ctls_size, ctl_flags are not valid after this call.
    UpdateNewRow(ctl_flags);

    // Add the column index
    if (full_column_indices_)
        da_put_u32(ctl_da_, x_start-1);
    else
        da_put_ul(ctl_da_, xs[0]);

    // Add deltas (if needed).
    if (xs_size > 1) {
        delta_bytes = DeltaSize_getBytes(delta_size);
        dst = dynarray_alloc_nr_aligned(ctl_da_, delta_bytes*(xs_size-1),
                                        delta_bytes);
        switch (delta_size) {
        case DELTA_U8:
            Copy((uint8_t  *) dst, &xs[1], xs_size-1);
            break;
        case DELTA_U16:
            Copy((uint16_t *) dst, &xs[1], xs_size-1);
            break;
        case DELTA_U32:
            Copy((uint32_t *) dst, &xs[1], xs_size-1);
            break;
        default:
            assert(false);
	    }
    }

    xs.clear();
    return;
}

template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
AddPattern(const Elem<IndexType, ValueType> &elem)
{
    uint8_t *ctl_flags, *ctl_size;
    long pat_size, pat_id;
    IndexType ucol;

    pat_size = elem.pattern->GetSize();
    pat_id = elem.pattern->GetPatternId();
    ctl_flags = (uint8_t *) dynarray_alloc_nr(ctl_da_, 2);
    *ctl_flags = GetFlag(pat_id, pat_size);
    ctl_size = ctl_flags + 1;
    assert(pat_size <= CTL_SIZE_MAX);
    *ctl_size = pat_size;
    UpdateNewRow(ctl_flags);

    if (full_column_indices_)
        ucol = elem.col;
    else
        ucol = elem.col - last_col_;
        
    LOG_DEBUG << "AddPattern ujmp " << ucol << "\n";

    if (full_column_indices_)
        da_put_u32(ctl_da_, ucol-1);
    else
        da_put_ul(ctl_da_, ucol);

    last_col_ = elem.pattern->ColIncreaseJmp(spm_->GetType(), elem.col);
    LOG_DEBUG << "last_col: " << last_col_ << "\n";
}

// return ujmp
template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
PreparePat(std::vector<IndexType> &xs, const Elem<IndexType, ValueType> &elem)
{
    if (xs.size() != 0)
        AddXs(xs);
}

} // end of csx namespace

#endif  // __cplusplus

#endif  // CSX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
