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
#include "spm_mt.h"

#ifdef __cplusplus
}
#endif

///< CSX matrix format
typedef struct {
    uint64_t nnz, ncols, nrows, ctl_size, row_start;
    double *values;
    uint8_t *ctl;
    long id_map[CTL_PATTERNS_MAX];
} csx_double_t;

typedef struct {
    csx_double_t *lower_matrix;
    double *dvalues;
} csx_double_sym_t;

///> Deletes all the fields of CSX format.
void DestroyCsx(csx_double_t *csx);
    
///> Deletes all the fields of CSX-Sym format.
void DestroyCsxSym(csx_double_sym_t *csx_sym);

// void DoPrint(void *csx);
// void PrintCsx(void *csx);

// #ifdef __cplusplus

// template<typename IndexType, typename ValueType>
// ValueType GetValueCsx(void *csx, IndexType row, IndexType col);
// template<typename IndexType, typename ValueType>
// bool SetValueCsx(void *csx, IndexType row, IndexType col, ValueType value);
// template<typename IndexType, typename ValueType>
// bool SearchValue(void *csx, IndexType row, IndexType col, ValueType& value, bool mode);

// #endif

#include "dynarray.h"

#define PID_DELTA_BASE 0

#ifdef __cplusplus  // This file is also included from csx_llvm_tmpl.c.

#include "SparsePartition.h"
#include "delta.h"

#ifdef SPM_NUMA
#   include <numa.h>
#   include "numa_util.h"
#   define DYNARRAY_CREATE  dynarray_create_numa
#else
#   define DYNARRAY_CREATE  dynarray_create
#endif

static bool debug = false;

namespace csx {

// template<typename IndexType, typename ValueType>
// struct csx_t {
//     IndexType nnz, ncols, nrows, ctl_size, row_start;
//     ValueType *values;
//     uint8_t *ctl;
//     long id_map[CTL_PATTERNS_MAX];
// };

// typedef csx_t<uint64_t, double> csx_double_t;

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
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
          ctl_da_(NULL),
          last_col_(0), empty_rows_(0) {}

#if 0   // SYM    
    CsxManager(SparsePartitionSym<IndexType, ValueType> *spm_sym)
        : spm_sym_(spm_sym),
          flag_avail_(0),
          row_jmps_(false),
          full_column_indices_(false),
          ctl_da_(NULL),
          last_col_(0), empty_rows_(0) {}
#endif  // SYM    

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
    csx_double_t *MakeCsx(bool symmetric);
    
    /**
     *  Transform the matrix owned by this manager into CSX form.
     *
     *  @return a handle to the newly created CSX-Sym matrix.
     */
#if 0   // SYM    
    csx_double_sym_t *MakeCsxSym();
#endif  // SYM    

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
               const Elem<IndexType, ValueType> *rend);
#if 0   // SYM    
    void DoSymRow(const Elem<IndexType, ValueType> *rstart,
                  const Elem<IndexType, ValueType> *rend);
#endif  // SYM    

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

    SparsePartition<IndexType, ValueType> *spm_;
    //SparsePartitionSym<IndexType, ValueType> *spm_sym_;
    uint8_t flag_avail_;     // Current available flags for pattern id mapping.
    bool row_jmps_;          // Whether or not row jumps included.
    bool full_column_indices_;  // use full 32-bit indices instead of deltas.
                                
    ValueType *values_;
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

template<typename T>
void Copy(T *dst, uint64_t *src, long nr_items)
{
    for (long i = 0; i < nr_items; i++){
        dst[i] = static_cast<T>(src[i]);
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
        typename CsxManager<IndexType, ValueType>::PatInfo pat_info(ret, 1, nnz);
        this->patterns[pattern_id] = pat_info;
    } else {
        ret = pi->second.flag;
        pi->second.npatterns++;
        pi->second.nr += nnz;
    }

    return ret;
}

#if 0   // SYM
csx_double_sym_t *CsxManager::MakeCsxSym()
{
    csx_double_sym_t *csx;
    double *diagonal = spm_sym_->GetDiagonal();
    uint64_t diagonal_size = spm_sym_->GetDiagonalSize();
    
    spm_ = spm_sym_->GetLowerMatrix();

#ifdef SPM_NUMA
    int cpu = sched_getcpu();
    if (cpu < 0) {
        perror("sched_getcpu() failed");
        exit(1);
    }

    int node = numa_node_of_cpu(cpu);
    if (node < 0) {
        perror("numa_node_of_cpu() failed");
        exit(1);
    }

    csx = (csx_double_sym_t *) numa_alloc_onnode(sizeof(csx_double_sym_t),
                                                 node);
    csx->dvalues = (double *) numa_alloc_onnode(diagonal_size * sizeof(double),
                                                node);
#else  
    csx = (csx_double_sym_t *) malloc(sizeof(csx_double_sym_t));
    csx->dvalues = (double *) malloc(diagonal_size * sizeof(double));
#endif

    for (uint64_t i = 0; i < diagonal_size; i++)
        csx->dvalues[i] = diagonal[i];
    csx->lower_matrix = MakeCsx(true);
    
    return csx;
}
#endif  // SYM    

template<typename IndexType, typename ValueType>
csx_double_t *CsxManager<IndexType, ValueType>::MakeCsx(bool symmetric)
{
    csx_double_t *csx;

#ifdef SPM_NUMA
    int cpu = sched_getcpu();
    if (cpu < 0) {
        perror("sched_getcpu() failed");
        exit(1);
    }

    int node = numa_node_of_cpu(cpu);
    if (node < 0) {
        perror("numa_node_of_cpu() failed");
        exit(1);
    }

    csx = (csx_double_t *) alloc_onnode(sizeof(csx_double_t), node);
    values_ = (ValueType *) alloc_onnode(sizeof(ValueType)*spm_->nr_nzeros_,
                                      node);
#else    
    csx = (csx_double_t *) malloc(sizeof(csx_double_t));
    values_ = (ValueType *) malloc(sizeof(ValueType)*spm_->nr_nzeros_);
#endif  // SPM_NUMA
    if (!csx || !values_) {
        std::cerr << __FUNCTION__ << ": malloc failed\n";
        exit(1);
    }

    // Be greedy with the initial capacity (equal to CSR col_ind size) 
    // to avoid realloc()'s.
    ctl_da_ = DYNARRAY_CREATE(sizeof(uint8_t), 512, 6*spm_->nr_nzeros_);
    csx->nnz = spm_->nr_nzeros_;
    csx->nrows = spm_->nr_rows_;
    csx->ncols = spm_->nr_cols_;
    csx->row_start = spm_->row_start_;
    values_idx_ = 0;
    new_row_ = false;		        // Do not mark first row.

    if (!symmetric) {
        for (IndexType i = 0; i < spm_->GetNrRows(); i++) {
            const Elem<IndexType, ValueType> *rbegin, *rend;
    
            rbegin = spm_->RowBegin(i);
            rend = spm_->RowEnd(i);
            if (debug)
                std::cerr << "MakeCsx(): row: " << i << "\n";
    
            if (rbegin == rend){		// Check if row is empty.
                if (debug)
                    std::cerr << "MakeCsx(): row is empty" << std::endl;
    
                if (new_row_ == false)
                    new_row_ = true;	// In case the first row is empty.
                else
                    empty_rows_++;
    
                continue;
            }
    
            DoRow(rbegin, rend);
            new_row_ = true;
        }
    } 
#if 0   // SYM    
    else {
        for (IndexType i = 0; i < spm_->GetNrRows(); i++) {
            const Elem<IndexType, ValueType> *rbegin, *rend;
    
            rbegin = spm_->RowBegin(i);
            rend = spm_->RowEnd(i);
            if (debug)
                std::cerr << "MakeCsx(): row: " << i << "\n";
    
            if (rbegin == rend){		// Check if row is empty.
                if (debug)
                    std::cerr << "MakeCsx(): row is empty" << std::endl;
    
                if (new_row_ == false)
                    new_row_ = true;	// In case the first row is empty.
                else
                    empty_rows_++;
    
                continue;
            }
    
            DoSymRow(rbegin, rend);
            new_row_ = true;
        }
    }
#endif  // SYM    

    csx->ctl_size = dynarray_size(ctl_da_);
    csx->ctl = (uint8_t *) dynarray_destroy(ctl_da_);
    ctl_da_ = NULL;
    assert(values_idx_ == spm_->nr_nzeros_);
    csx->values = values_;
    values_ = NULL;
    values_idx_ = 0;

    AddMappings(csx->id_map);
    // for (int i=0; csx->id_map[i] != -1; i++)
    //     std::cout << csx->id_map[i] << " ";
    // std:: cout << std::endl;
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

/*
 *  Ctl Rules
 *  1. Each unit leaves the x index at the last element it calculated on the
 *     current row.
 *  2. Size is the number of elements that will be calculated.
 */
template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
DoRow(const Elem<IndexType, ValueType> *rbegin,
      const Elem<IndexType, ValueType> *rend)
{
    std::vector<IndexType> xs;

    last_col_ = 1;
    for (const Elem<IndexType, ValueType> *spm_elem = rbegin;
         spm_elem < rend; spm_elem++) {
        if (debug)
            std::cerr << "\t" << *spm_elem << "\n";

        // Check if this element contains a pattern.
        if (spm_elem->pattern != NULL) {
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

#if 0   // SYM
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

    last_col_ = 1;
    for ( ; spm_elem < rend && spm_elem->col < spm_->GetRowStart() + 1; 
         spm_elem++) {
        if (debug)
            std::cerr << "\t" << *spm_elem << "\n";

        // Check if this element contains a pattern.
        if (spm_elem->pattern != NULL) {
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
        if (debug)
            std::cerr << "\t" << *spm_elem << "\n";

        // Check if this element contains a pattern.
        if (spm_elem->pattern != NULL) {
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
#endif  // SYM    

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
        
    if (debug)
        std::cerr << "AddPattern ujmp " << ucol << "\n";

    if (full_column_indices_)
        da_put_u32(ctl_da_, ucol-1);
    else
        da_put_ul(ctl_da_, ucol);

    last_col_ = elem.pattern->ColIncreaseJmp(spm_->GetType(), elem.col);
    if (debug)
        std::cerr << "last_col:" << last_col_ << "\n";
}

// return ujmp
template<typename IndexType, typename ValueType>
void CsxManager<IndexType, ValueType>::
PreparePat(std::vector<IndexType> &xs, const Elem<IndexType, ValueType> &elem)
{
    if (xs.size() != 0)
        AddXs(xs);
}

}   // end of csx namespace

/* Get/Set utilities of CSX */

/*
 * Returns the value of the corresponding element (0 if element doesn't exist)
 * Indexes row/col are assumed to be 1-based.
 */
template<typename IndexType, typename ValueType>
ValueType GetValueCsx(void *spm, IndexType row, IndexType col)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
    ValueType value = 0;

    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        SearchValue(spm_mt->spm_threads[i].spm, row, col, value, false);
    }
    return value;
}

/*
 * Sets the value of the corresponding element.
 * Returns false if the element doesn't exist.
 * Indexes row/col are assumed to be 1-based.
 */
template<typename IndexType, typename ValueType>
bool SetValueCsx(void *spm, IndexType row, IndexType col, ValueType value)
{
	spm_mt_t *spm_mt = (spm_mt_t *) spm;
    bool found = false;

    for (size_t i = 0; i < spm_mt->nr_threads; i++) {
        found = SearchValue(spm_mt->spm_threads[i].spm, row, col, value, true);
    }
    return found;
}

#define ALIGN(buf,a) (uint8_t *) (((unsigned long) (buf) + (a-1)) & ~(a-1))
void align_ptr(uint8_t **ctl, int align)
{
    *ctl = ALIGN(*ctl, align);
}

/*
 * Searches for the corresponding element.
 * Returns false if the element doesn't exist.
 * Indexes row/col are assumed to be 1-based.
 */
template<typename IndexType, typename ValueType>
bool SearchValue(void *spm, IndexType row, IndexType col, ValueType &value,
                 bool mode)
{
    csx_double_t *csx = (csx_double_t *) spm;
	ValueType *v = csx->values;
	uint8_t *ctl = csx->ctl;
	uint8_t *ctl_end = ctl + csx->ctl_size;
	uint8_t size, flags;
	uint8_t patt_id, type, delta;
    IndexType cur_row = csx->row_start;
    IndexType ucol = 0;
    uint8_t row_dim, col_dim;
    uint64_t i_end;

	do {
        flags = *ctl++;
		size = *ctl++;
		if (test_bit(&flags, CTL_NR_BIT)) {
            if (test_bit(&flags, CTL_RJMP_BIT)) {
                cur_row += ul_get(&ctl);
            } else {
                cur_row++;
            }
            ucol = 0;
		}
        if (false)  //(full_column_indices)
            ucol = u32_get(&ctl);// + last_col;
        else
            ucol += ul_get(&ctl);

		patt_id = flags & CTL_PATTERN_MASK;
        type = csx->id_map[patt_id] / CSX_PID_OFFSET;
        delta = csx->id_map[patt_id] % CSX_PID_OFFSET;  //other_dim in case it is a block

//        if (type > BLOCK_TYPE_START && type < BLOCK_COL_START) {
        if (type >= 6 && type <= 13) {
            row_dim = type - BLOCK_TYPE_START;
            col_dim = delta;
//        } else if (type > BLOCK_COL_START && type <= BLOCK_TYPE_END) {
        } else if (type >= 15 && type <= 22) {
            col_dim = type - BLOCK_COL_START;
            row_dim = delta;
        }

        // std::cout << "row:" << cur_row << " "; 
        // std::cout << "col:" << ucol << " ";
        // printf("type:%"PRIu8 " ", type);
        // printf("delta:%"PRIu8 " ", delta);
        // printf("size:%"PRIu8 " ", size);
        // std::cout << std::endl;

        if ((row - 1) == cur_row) {
            if ((col - 1) == ucol) {
                if (mode)
                    *v = value;
                else
                    value = *v;
                return true;
            } else if ((col - 1) < ucol) {
                return false;
            } else {
                switch (type) {
                case 0:  // DELTAS
                    if (delta == 8) {
                        for (uint8_t i = 1; i < size; i++) {
                            v++;
                            ucol += u8_get(&ctl);
                            if ((col - 1) == ucol) {
                                if (mode)
                                    *v = value;
                                else
                                    value = *v;
                                return true;
                            }
                        }
                    } else if (delta == 16) {
                        align_ptr(&ctl, 2);
                        for (uint8_t i = 1; i < size; i++) {
                            v++;
                            ucol += u16_get(&ctl);
                            if ((col - 1) == ucol) {
                                if (mode)
                                    *v = value;
                                else
                                    value = *v;
                                return true;
                            }
                        }
                    } else if (delta == 32){
                        align_ptr(&ctl, 4);
                        for (uint8_t i = 1; i < size; i++) {
                            v++;
                            ucol += u32_get(&ctl);
                            if ((col - 1) == ucol) {
                                if (mode)
                                    *v = value;
                                else
                                    value = *v;
                                return true;
                            }
                        }
                    } else {
                        align_ptr(&ctl, 8);
                        for (uint8_t i = 1; i < size; i++) {
                            v++;
                            ucol += u64_get(&ctl);
                            if ((col - 1) == ucol) {
                                if (mode)
                                    *v = value;
                                else
                                    value = *v;
                                return true;
                            }
                        }
                    }
                    v++;
                    break;
                case 1: //  HORIZONTAL
                    i_end = size * delta;
                    if ((col - 1) < (ucol + size * delta)) {
                        for (uint64_t i = delta; i < i_end; i += delta) {
                            v++;
                            if ((col - 1) == (ucol + i)) {
                                if (mode)
                                    *v = value;
                                else
                                    value = *v;
                                return true;
                            }
                        }
                        v++;
                    } else {
                        v += size;
                    }
                    ucol += i_end - delta;
                    break;
                case 2 ... 4:   //VERTICAL, DIAGONAL, REV_DIAGONAL: can't possibly be in here
                    v += size;
                    break;
                case 6 ... 13:  //BLOCK_R1 ... BLOCK_COL_START - 1  
                    if ((col - 1) < (ucol + col_dim)) {
                        v += row_dim * (col - 1 - ucol);
                        if (mode)
                            *v = value;
                        else
                            value = *v;
                        return true;
                    } else {
                        v += size;
                    }
                    break;
                case 15 ... 22:  //BLOCK_COL_START ... BLOCK_TYPE_END - 1
                    if ((col - 1) < (ucol + col_dim)) {
                        v += (col-1) - ucol;
                        if (mode)
                            *v = value;
                        else
                            value = *v;
                        return true;
                    } else {
                        v += size;
                    }
                    break;
                default:
                    assert(0 && "unknown pattern type");
                }
            }
        } else if ((row - 1) < cur_row) {
            return false;
        } else {
            switch (type) {
            case 0: //NONE
                if (delta == 8) {
                    for (uint8_t i = 1; i < size; i++) {
                        ucol += u8_get(&ctl);
                        v++;
                    }
                } else if (delta == 16) {
                    align_ptr(&ctl, 2);
                    for (uint8_t i = 1; i < size; i++) {
                        ucol += u16_get(&ctl);
                        v++;
                    }
                } else if (delta == 32){
                    align_ptr(&ctl, 4);
                    for (uint8_t i = 1; i < size; i++) {
                        ucol += u32_get(&ctl);
                        v++;
                    }
                } else {
                    align_ptr(&ctl, 8);
                    for (uint8_t i = 1; i < size; i++) {
                        ucol += u64_get(&ctl);
                        v++;
                    }
                }
                v++;
                break;
            case 1: //HORIZONTAL: can't possibly be in here
                v += size;
                ucol += (size - 1) * delta;
                break;
            case 2: //VERTICAL
                if (((col - 1) == ucol) &&
                    ((row - 1) < (cur_row + size * delta))) {
                    i_end = size * delta;
                    for (uint64_t i = delta; i < i_end; i += delta) {
                        v++;
                        if ((row - 1) == (cur_row + i)) {
                            if (mode)
                                *v = value;
                            else
                                value = *v;
                            return true;
                        }
                    }
                    v++;
                } else {                    
                    v += size;
                }
                break;
            case 3: //DIAGONAL
                if ((col - 1) <= ucol) {  // can't possibly be in here
                    v += size;
                } else if ((col - 1) <= (ucol + (size - 1) * delta)) {
                    i_end = size * delta;
                    for (uint64_t i = delta; i < i_end; i += delta) {
                        v++;
                        if (((row - 1) == (cur_row + i)) && ((col - 1) == (ucol + i))) {
                            if (mode)
                                *v = value;
                            else
                                value = *v;
                            return true;
                        }                        
                    }
                    v++;
                } else {
                    v += size;
                } 
                break;
            case 4: //REV_DIAGONAL
                if ((col - 1) >= ucol) {  // can't possibly be in here
                    v += size;
                } else if ((col - 1) >= (ucol - (size - 1) * delta)) {
                    i_end = size * delta;
                    for (uint64_t i = delta; i < i_end; i += delta) {
                        v++;
                        if (((row - 1) == (cur_row + i)) && ((col - 1) == (ucol - i))) {
                            if (mode)
                                *v = value;
                            else
                                value = *v;
                            return true;
                        }                        
                    }
                    v++;
                } else {
                    v += size;
                } 
                break;
            case 6 ... 13:  //BLOCK_R1 ... BLOCK_COL_START - 1
                if (((col - 1) < (ucol + col_dim))
                    && ((row -1) < (cur_row + row_dim))) {
                    for (uint64_t i = 0; i < col_dim; i++) {
                        for (uint64_t j = 0; j < row_dim; j++) {
                            if (((row - 1) == (cur_row + j)) && ((col - 1) == (ucol + i))) {
                                if (mode)
                                    *v = value;
                                else
                                    value = *v;
                                return true;
                            }                        
                            v++;
                        }
                    } 
                } else {                    
                    v += size;
                }
                break;
            case 15 ... 22: //BLOCK_COL_START ... BLOCK_TYPE_END - 1
                if (((col - 1) < (ucol + col_dim))
                    && ((row -1) < (cur_row + row_dim))) {
                    for (uint64_t i = 1; i < row_dim; i++) {
                        v += col_dim;
                        for (uint64_t j = 0; j < col_dim; j++) {
                            if (((row - 1) == (cur_row + i)) && ((col - 1) == (ucol + j))) {
                                if (mode)
                                    *(v + j) = value;
                                else
                                    value = *(v + j);
                                return true;
                            }                        
                        }
                    } 
                    v += col_dim;
                } else {                    
                    v += size;
                }
                break;
            default:
                assert(0 && "unknown pattern type");
            }
        }
	} while (ctl < ctl_end);

    return false;
}

#endif  // __cplusplus

#endif  // CSX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
                    // if (size > 1) {
                    //     size_t count = 1;
                    //     if (delta == 16) {
                    //         align_ptr(&ctl, 2);
                    //     } else if (delta == 32) {
                    //         align_ptr(&ctl, 4);
                    //     } else if (delta == 64){
                    //         align_ptr(&ctl, 8);
                    //     }
                    //     do {
                    //         v++;
                    //         count++;
                    //         if (delta == 8) {
                    //             ucol += u8_get(&ctl);
                    //         } else if (delta == 16) {
                    //             ucol += u16_get(&ctl);
                    //         } else if (delta == 32) {
                    //             ucol += u32_get(&ctl);
                    //         } else {
                    //             ucol += u64_get(&ctl);
                    //         }
                    //         if ((col - 1) == tmp_col) {
                    //             if (mode)
                    //                 *v = value;
                    //             else
                    //                 value = *v;
                    //             return true;
                    //         }
                    //         if ((col - 1) == tmp_col) {
                    //             if (mode)
                    //                 *v = value;
                    //             else
                    //                 value = *v;
                    //             return true;
                    //         }
                    //     } while (count < size);
                    //     v++;
                    // } else {
                    //     last_col = ucol;
                    //     v++;
                    // }

                // if (size > 1) {
                //     IndexType tmp_col = ucol;
                //     size_t count = 1;
                //     do {
                //         v++;
                //         count++;
                //         if (delta8) {
                //             delta = u8_get(&ctl);
                //             tmp_col += delta;
                //         }
                //         else if (delta16) {
                //             delta = u8_get(&ctl);
                //             tmp_col += delta;
                //         }
                //         else {
                //             delta = u32_get(&ctl);
                //             tmp_col += delta;
                //         }
                //     } while (count < size);
                //     last_col = tmp_col;
                //     v++;
                // } else {
                //     last_col = ucol;
                //     v++;
                // }
                    // if ((col - 1) < (ucol + size * delta)) {
                    //     IndexType tmp_col = ucol;
                    //     size_t count = 1;
                    //     do {
                    //         tmp_col += delta;
                    //         v++;
                    //         if ((col - 1) == tmp_col) {
                    //             if (mode)
                    //                 *v = value;
                    //             else
                    //                 value = *v;
                    //             return true;
                    //         }
                    //         count++;
                    //     } while (count < size);
                    // } else {
                    //     v += size;
                    // }

                    // IndexType tmp_row = cur_row;
                    // size_t count = 1;
                    // do {
                    //     tmp_row += delta;
                    //     v++;
                    //     if ((row - 1) == tmp_row) {
                    //         if (mode)
                    //             *v = value;
                    //         else
                    //             value = *v;
                    //         return true;
                    //     }   
                    //     count++;
                    // } while (count < size);
                    // v++;


                    // IndexType tmp_row = cur_row;
                    // IndexType tmp_col = ucol;
                    // size_t count = 1;
                    // do {
                    //     tmp_row += delta;
                    //     tmp_col += delta;
                    //     v++;
                    //     if (((col-1) == tmp_col) && ((row - 1) == tmp_row)) {
                    //         if (mode)
                    //             *v = value;
                    //         else
                    //             value = *v;
                    //         return true;
                    //     }
                    //     count++;
                    // } while (count < size);
                    // v++;


                    // IndexType tmp_row = cur_row;
                    // IndexType tmp_col = ucol;
                    // size_t count = 1;
                    // do {
                    //     tmp_col -= delta;
                    //     tmp_row += delta;
                    //     v++;
                    //     if (((col-1) == tmp_col) && ((row - 1) == tmp_row)) {
                    //         if (mode)
                    //             *v = value;
                    //         else
                    //             value = *v;
                    //         return true;
                    //     }
                    //     count++;
                    // } while (count < size);
                    // v++;



                    // v += delta;
                    // IndexType tmp_row = cur_row + 1;
                    // IndexType tmp_col = ucol;
                    // size_t count = delta;
                    // do {
                    //     if (((col-1) == tmp_col) && ((row - 1) == tmp_row)) {
                    //         if (mode)
                    //             *v = value;
                    //         else
                    //             value = *v;
                    //         return true;
                    //     }
                    //     tmp_col++;
                    //     if (tmp_col == (ucol + delta)) {
                    //         tmp_col = ucol;
                    //         tmp_row++;
                    //     }
                    //     v++;
                    //     count++;
                    // } while (count < size);


                    // v += first_dim;
                    // IndexType tmp_row = cur_row + 1;
                    // IndexType tmp_col = ucol;
                    // size_t count = first_dim;
                    // do {
                    //     if (((col-1) == tmp_col) && ((row - 1) == tmp_row)) {
                    //         if (mode)
                    //             *v = value;
                    //         else
                    //             value = *v;
                    //         return true;
                    //     }
                    //     tmp_col++;
                    //     if (tmp_col == (ucol + first_dim)) {
                    //         tmp_col = ucol;
                    //         tmp_row++;
                    //     }
                    //     v++;
                    //     count++;
                    // } while (count < size);
