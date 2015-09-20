/*
 * Copyright (C) 2009-2015, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2015, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file SparsePartition.hpp
 * \brief Internal representation of sparse matrice partition
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2015
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_SPARSE_PARTITION_HPP
#define SPARSEX_INTERNALS_SPARSE_PARTITION_HPP

#include <sparsex/internals/Allocators.hpp>
#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/DynamicArray.hpp>
#include <sparsex/internals/Element.hpp>
#include <sparsex/internals/Utility.hpp>
#include <sparsex/internals/Xform.hpp>
#include <cassert>
#include <vector>
#include <iterator>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace sparsex::utilities;

namespace sparsex {
namespace csx {

template<typename IndexType, typename ValueType>
class SparsePartitionSym;

/**
 *  Internal representation of a sparse matrix. This class can also represent a
 *  sub-matrix of a larger one, starting at a specific row of the original
 *  matrix.
 */
template<typename IndexType, typename ValueType>
class SparsePartition
{
public:
    typedef IndexType idx_t;
    typedef ValueType val_t;
#if SPX_USE_NUMA
    typedef NumaAllocator MemoryAllocatorImpl;
#else
    typedef StdAllocator MemoryAllocatorImpl;
#endif

    SparsePartition() :
        nr_rows_(0),
        nr_cols_(0),
        nr_nzeros_(0),
        type_(Encoding::None),
        elems_size_(0),
        rowptr_(NULL),
        rowptr_size_(0),
        elems_mapped_(false),
        alloc_(MemoryAllocatorImpl::GetInstance())
    {}

    ~SparsePartition()
    {
        alloc_.Destroy(rowptr_, rowptr_size_);
        if (!elems_mapped_)
            elems_.clear();
    };

    size_t GetNrRows() const
    {
        return nr_rows_;
    }

    size_t GetNrCols() const
    {
        return nr_cols_;
    }

    size_t GetNrNonzeros() const
    {
        return nr_nzeros_;
    }

    size_t GetNrDeltas() const
    {
        return nr_deltas_;
    }

    Encoding::Type GetType() const
    {
        return type_;
    }

    IndexType GetRow(IndexType idx) const
    {
        return rowptr_[idx];
    }
    
    IndexType GetRowStart() const
    {
        return row_start_;
    }

    size_t GetRowptrSize() const
    {
        return rowptr_size_;
    }

    vector<Element<IndexType, ValueType> > &GetElems()
    {
        return elems_;
    }

    void SetNrRows(size_t nr_rows)
    {
        nr_rows_ = nr_rows;
    }

    void SetNrCols(size_t nr_cols)
    {
        nr_cols_ = nr_cols;
    }

    void SetNrNonzeros(size_t nr_nzeros)
    {
        nr_nzeros_ = nr_nzeros;
    }

    void SetType(Encoding::Type type)
    {
        type_ = type;
    }
    
    void SetRowStart(IndexType row_start)
    {
        row_start_ = row_start;
    }
    
    void SetElemsSize(size_t elems_size)
    {
        elems_size_ = elems_size;
    }

    void SetElem(IndexType idx, Element<IndexType, ValueType> &&elem)
    {
        elems_[idx] = elem;
    }

    void AddNrDeltas(size_t nr)
    {
        nr_deltas_ += nr;
    }

    Element<IndexType, ValueType>& operator[](const size_t idx)
    {
        assert(idx >= 0 && idx < nr_nzeros_ && "Index out of bounds!");
        return elems_[idx];
    }

    class Builder;

    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start point iterators start.
     *  @param pnts_end   point iterators end.
     *  @param row_start  first row of the matrix.
     *  @param limit      limit for the number of elements.
     *  @param nr_elems   estimation for the number of elements.
     *  @param nr_rows    size for rows allocation.
     *  @return           number of elements actually allocated.
     */
    template <typename IterT>
    IndexType SetElems(IterT &pnts_start, const IterT &pnts_end,
                       IndexType row_start, size_t limit,
                       size_t nr_elems, size_t nr_rows);

    /**
     *  Function for constructing the new rowptr array using point
     *  iterators
     *
     *  @param pnts_start point iterators start.
     *  @param pnts_end   point iterators end.
     *  @param row_start  first row of the matrix.
     *  @param nr_rows    size for rows allocation.
     */
    template <typename IterT>
    void SetRowptr(IterT &pnts_start, const IterT &pnts_end,
                   size_t nr_rows);

    /**
     *  Printing methods.
     */
    void Print(ostream &out = cout);
    void PrintRows(ostream &out = cout);
    void PrintStats(ostream &out = cout);

    template<typename I, typename V>
    friend ostream &operator<<(ostream &, const SparsePartition<I, V> &);

    IndexType FindNewRowptrSize(Encoding::Type t);

    /**
     *  Transforms matrix to specific type.
     *
     *  @param type type of pattern to transform to.
     *  @param rs   line to start the transformation. If it is zero, it is the
     *              the beggining of matrix.
     *  @param re   line to end the transformation. If it is zero, it is the
     *              the end of matrix.
     */
    void Transform(Encoding::Type type, IndexType rs = 0, IndexType re = 0);

    /**
     *  Extracts a window (sub-matrix) from this matrix. This methods copies the
     *  nonzero elements of the original matrix to the window.
     *
     *  @param rs     the row in this matrix, from where the window will be
     *                extracted.
     *  @param length number of rows to include in the window.
     *  @return       the selected window.
     *  @see GetWindow()
     */
    SparsePartition *ExtractWindow(IndexType rs, IndexType length);

    /**
     *  Extracts a window (sub-matrix) from this matrix using a backing
     *  store. This method uses the nonzero elements of this matrix without
     *  performing any copy.
     *
     *  @param rs     the row in this matrix, from where the window will be
     *                extracted.
     *  @param length number of rows to include in the window.
     *  @return       the selected window.
     *  @see ExtractWindow()
     */
    SparsePartition *GetWindow(IndexType rs, IndexType length);

    /**
     *  Puts a (possibly processed) window of the matrix back to this matrix.
     *
     *  @param window the (possibly processed) window.
     */
    void PutWindow(SparsePartition *window);
 
    // SparsePartition's iterator
    class iterator;
    iterator begin(IndexType ridx)
    {
        return iterator(this, ridx);
    }

    iterator end(IndexType ridx)
    {
        struct iterator::end_tag end;
        return iterator(this, ridx, end);
    }

    friend class SparsePartitionSym<IndexType, ValueType>;

private:
    size_t nr_rows_, nr_cols_, nr_nzeros_;
    Encoding::Type type_;
    vector<Element<IndexType, ValueType> > elems_;
    size_t elems_size_;
    IndexType *rowptr_;
    IndexType rowptr_size_;
    IndexType row_start_;    /* Row of the original matrix, where this
                                sub-matrix starts. */
    bool elems_mapped_;
    size_t nr_deltas_;
    MemoryAllocator &alloc_;
};

/**
 *  Utility class for building a SparsePartition.
 */
template<typename IndexType, typename ValueType>
class SparsePartition<IndexType, ValueType>::Builder
{
public:
#if SPX_USE_NUMA
    typedef DynamicArray<IndexType,
                         reallocator<IndexType,
                                     NumaAllocator> > DynamicIndexArray;
#else
    typedef DynamicArray<IndexType> DynamicIndexArray;
#endif

    /**
     *  Constructs a Builder object.
     *
     *  @param sp       the matrix to be built.
     *  @param nr_rows  an initial guess of the number of rows.
     */
    Builder(SparsePartition *sp, size_t nr_rows = 0);

    ~Builder();

    /**
     *  Allocate elements for the next row.
     *
     *  @param rdiff number of empty rows between the current and the new row.
     *  @param value number of elements for the current row.
     */
    void NewRow(IndexType rdiff = 1, IndexType value = 0);

    /**
     *  Finalizes the construction of the new matrix.
     */
    void Finalize(size_t size);

public:
    SparsePartition *sp_;
    DynamicIndexArray da_rowptr_;
};

/**
 *  Iterators for moving between elements inside the matrix.
 */
template<typename IndexType, typename ValueType>
class SparsePartition<IndexType, ValueType>::iterator
    : public std::iterator<std::forward_iterator_tag,
                           Element<IndexType, ValueType> >
{
public:
    struct end_tag {};
    iterator(SparsePartition *sp, IndexType row_idx);
    iterator(SparsePartition *sp, IndexType row_idx, struct end_tag);
    bool operator==(const iterator &pi);
    bool operator!=(const iterator &pi);
    void operator++();
    Element<IndexType, ValueType> &operator*();

private:
    SparsePartition *sp_;
    IndexType row_idx_;
    IndexType elem_idx_;
};

/**
 *  Internal representation of a sparse symmetric matrix. This class can also 
 *  represent a sub-matrix of a larger one.
 */
template<typename IndexType, typename ValueType>
class SparsePartitionSym
{
private:
    SparsePartition<IndexType, ValueType> *lower_matrix_; // lower triangular
    vector<ValueType> diagonal_;   
    IndexType diagonal_size_;
    SparsePartition<IndexType, ValueType> *m1_; // Matrix that contains the
                                                // elems of SpmSym for which
                                                // column of the element is
                                                // smaller than the row start.

    SparsePartition<IndexType, ValueType> *m2_; // Matrix that contains the rest
                                                // of the elements.
    MemoryAllocator &alloc_;
#if SPX_USE_NUMA
    typedef NumaAllocator MemoryAllocatorImpl;
#else
    typedef StdAllocator MemoryAllocatorImpl;
#endif
    
public:
    typedef IndexType idx_t;
    typedef ValueType val_t;

    SparsePartitionSym() 
        : lower_matrix_(NULL),
          m1_(NULL),
          m2_(NULL),
          alloc_(MemoryAllocatorImpl::GetInstance())
    {
        lower_matrix_ = new SparsePartition<IndexType, ValueType>;
        m1_ = new SparsePartition<IndexType, ValueType>;
        m2_ = new SparsePartition<IndexType, ValueType>;
    }

    ~SparsePartitionSym()
    {
        delete lower_matrix_;
        diagonal_.clear();
    }
    
    SparsePartition<IndexType, ValueType> *GetLowerMatrix()
    {
        return lower_matrix_;
    }
    
    SparsePartition<IndexType, ValueType> *GetFirstMatrix()
    {
        return m1_;
    }
    
    SparsePartition<IndexType, ValueType> *GetSecondMatrix()
    {
        return m2_;
    }
    
    vector<ValueType> *GetDiagonal()
    {
        return &diagonal_;
    }
    
    IndexType GetNrRows()
    {
        return diagonal_size_;
    }

    IndexType GetDiagonalSize()
    {
        return diagonal_size_;
    }
    
    void SetDiagonalSize(IndexType size)
    {
        diagonal_size_ = size;
    }
    
    void SetDiagonal(ValueType *diagonal)
    {
        diagonal_ = diagonal;
    }

    void SetNrRows(size_t nr_rows)
    {
        lower_matrix_->nr_rows_ = nr_rows;
    }

    void SetNrCols(size_t nr_cols)
    {
        lower_matrix_->nr_cols_ = nr_cols;
    }

    void SetNrNonzeros(size_t nnz)
    {
        lower_matrix_->nr_nzeros_ = nnz - diagonal_size_;
    }

    void SetRowStart(IndexType row_start)
    {
        lower_matrix_->row_start_ = row_start;
    }
    
    void SetType(Encoding::Type type)
    {
        lower_matrix_->type_ = type;
    }

    size_t GetRowptrSize() const
    {
        return lower_matrix_->rowptr_size_;
    }

    /**
     *  Divide the matrix (m) to two submatrices (m1 and m2).
     */
    void DivideMatrix();

    /**
     *  Merge the two submatrices (m1 and m2) into one matrix (m).
     */
    void MergeMatrix();
    
    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start   point iterators start.
     *  @param pnts_end     point iterators end.
     *  @param row_start    first row of the matrix.
     *  @param limit        limit for the number of elements.
     *  @param nr_elems     size for elems allocation.
     *  @param nr_rows      size for rows allocation.
     *  @return             number of elements allocated.
     */
    template <typename IterT>
    IndexType SetElems(IterT &pnts_start, const IterT &pnts_end,
                       IndexType row_start, size_t limit,
                       size_t nr_elems, size_t nr_rows);

    void PrintDiagElems(ostream &out = cout);
};

/* Helper functions */

template<typename IndexType, typename ValueType>
void TestMMF(SparsePartition<IndexType, ValueType> *sp, const char *mmf_file);

/*
 *  SparsePartition class implementation
 */
template<typename IndexType, typename ValueType> template<typename IterT>
IndexType SparsePartition<IndexType, ValueType>::
SetElems(IterT &pi, const IterT &pnts_end, IndexType row_start,
         size_t limit, size_t nr_elems, size_t nr_rows)
{
    SparsePartition::Builder Bld(this, nr_rows);
    IndexType row_prev = 1;
    size_t elem_cnt = 0;

    elems_.reserve(nr_elems);
    // elems_.reserve(limit);
    for (; pi != pnts_end; ++pi) {
        Element<IndexType, ValueType> elem(*pi);
        IndexType row = elem.GetRow() - row_start + 1;
        if (row != row_prev) {
            assert(row > row_prev);
            if (limit && elem_cnt >= limit)
                break;	   
            Bld.NewRow(row - row_prev, elem_cnt);
            row_prev = row;
        }
    
        // Element's row must be set to the new value
        elem.Transform(make_pair(row, elem.GetCol()));
        elems_.push_back(move(elem));
        elem_cnt++;
    }

    // elems_'s size will never increase henceforth, so adjust capacity
    elems_.shrink_to_fit();
    elems_size_ = elems_.size();
    Bld.Finalize(elems_size_);
    return elems_size_;
}

template<typename IndexType, typename ValueType> template<typename IterT>
void SparsePartition<IndexType, ValueType>::
SetRowptr(IterT &pi, const IterT &pnts_end, size_t nr_rows)
{
    SparsePartition::Builder Bld(this, nr_rows);
    IndexType row_prev = 1;
    size_t elem_cnt = 0;

    for (; pi != pnts_end; ++pi) {
        IndexType row = (*pi).GetRow();
        if (row != row_prev) {
            assert(row > row_prev);
            Bld.NewRow(row - row_prev, elem_cnt);
            row_prev = row;
        }
      
        elem_cnt++;
    }

    Bld.Finalize(elem_cnt);
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::Print(std::ostream &out)
{
    iterator p_start = begin(0);
    iterator p_end = end(rowptr_size_ - 2);
    for (auto p = p_start; p != p_end; ++p)
        out << " " << (*p);

    out << endl;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::PrintRows(ostream &out) {
    cout << "Row Ptr: ";
    for (IndexType i = 0; i < rowptr_size_; i++)
        cout << rowptr_[i] << " ";
    cout << endl;
}

// FIXME
template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::PrintStats(ostream& out)
{
    IndexType nr_rows_with_patterns;
    IndexType nr_patterns, nr_patterns_before;
    IndexType nr_nzeros_block;
    IndexType nr_transitions;
    IndexType nr_xform_patterns[Encoding::Max];

    nr_rows_with_patterns = rowptr_size_ - 1;
    nr_patterns = 0;
    nr_nzeros_block = 0;
    nr_transitions = 0;
    memset(nr_xform_patterns, 0, sizeof(nr_xform_patterns));
    for (size_t i = 0; i < rowptr_size_ - 1; i++) {
        IndexType pt_size, pt_size_before;
        Encoding::Type pt_type, pt_type_before;

        nr_patterns_before = nr_patterns;
        SparsePartition<IndexType, ValueType>::iterator ei = begin(i);
        SparsePartition<IndexType, ValueType>::iterator ee = end(i);
        if (ei->IsPattern()) {
            pt_size_before = pt_size = ei->GetSize();
            pt_type_before = pt_type = ei->GetType();
        } else {
            pt_size_before = pt_size = 0;
            pt_type_before = pt_type = Encoding::None;
        }

        for (; ei != ee; ++ei) {
            if (ei->IsPattern()) {
                ++nr_patterns;
                pt_size = ei->GetSize();
                pt_type = ei->GetType();
                nr_nzeros_block += pt_size;
                if (pt_type != pt_type_before ||
                    (pt_size_before && pt_size != pt_size_before))
                    ++nr_transitions;
                    
                ++nr_xform_patterns[ei->pattern->GetType()];
                pt_size_before = pt_size;
            } else {
                pt_type = Encoding::None;
                ++nr_xform_patterns[Encoding::None];
                if (pt_type != pt_type_before)
                    ++nr_transitions;
            }

            pt_type_before = pt_type;
        }

        if (nr_patterns == nr_patterns_before)
            --nr_rows_with_patterns;
    }

    size_t nr_encoded_types = 0;
    for (Encoding::Type t = Encoding::Horizontal; t < Encoding::Max; ++t)
        if (nr_xform_patterns[t]) {
            ++nr_encoded_types;
            Encoding enc(t);
            out << enc << ": " << nr_xform_patterns[t] << "\n";
        }

    out << "Encoded types = " << nr_encoded_types << ", "
        << "Avg patterns/row = "
        << nr_patterns / (double) nr_rows_with_patterns << ", "
        << "Avg nonzeros/pattern = "
        << nr_nzeros_block / (double) nr_patterns << ", "
        << " Avg pattern transitions/row = "
        << nr_transitions / (double) nr_rows_with_patterns
        << "\n";
}

template<typename IndexType, typename ValueType>
IndexType SparsePartition<IndexType, ValueType>::
FindNewRowptrSize(Encoding::Type t)
{
    Encoding e(t);
    if (t == Encoding::Horizontal)
        return nr_rows_ + 1;
    else if (t == Encoding::Vertical)
        return nr_cols_ + 1;
    else if (t == Encoding::Diagonal || t == Encoding::AntiDiagonal)
        return nr_rows_ + nr_cols_ + 1;
    else if (e.IsBlockRow())
        return nr_rows_ / e.GetBlockAlignment() + 2;
    else if (e.IsBlockCol())
        return nr_cols_ / e.GetBlockAlignment() + 2;
    else
        return 0;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::Transform(
    Encoding::Type t, IndexType rs, IndexType re)
{
    if (type_ == t)
        return;

    if (!re)
        re = rowptr_size_ - 2;

    typename vector<Element<IndexType, ValueType> >::iterator e0, ee, es;
    typename TransformFn<IndexType>::type xform_fn =
        GetXformFn<IndexType>(type_, t);

    // Transform element coordinates
    e0 = elems_.begin();
    ee = e0 + elems_size_;
    for (; e0 != ee; ++e0) {
        (*e0).Transform(xform_fn(make_pair((*e0).GetRow()-rs, (*e0).GetCol()),
                                 nr_rows_, nr_cols_));
    }

    Encoding e(t);
    Encoding e_(type_);
    if (((type_ == Encoding::Horizontal || e_.IsBlockRow()) &&
         (t == Encoding::Horizontal || e.IsBlockRow())) ||
        ((type_ == Encoding::Vertical || e_.IsBlockCol()) &&
         (t == Encoding::Vertical || e.IsBlockCol()))) {
        IndexType old_block_align, new_block_align, k;

        if (type_ == Encoding::Horizontal || type_ == Encoding::Vertical)
            old_block_align = 1;
        else
            old_block_align = e_.GetBlockAlignment();

        if (t == Encoding::Horizontal || t == Encoding::Vertical)
            new_block_align = 1;
        else
            new_block_align = e.GetBlockAlignment();

        k = lcm(old_block_align, new_block_align);
        k /= old_block_align;
        re = rowptr_size_ - 1;
        
        e0 = elems_.begin();
        ee = e0;
        for (IndexType i = k; i < re; i += k) {
            es = ee;
            ee += rowptr_[i] - rowptr_[i-k];
            sort(es, ee);
        }

        es = ee;
        ee = elems_.begin() + elems_size_;
        sort(es, ee);
    } else {
        e0 = elems_.begin();
        ee = e0 + elems_size_;
        sort(e0, ee);
    }

    if (elems_size_)
        SetRowptr(e0, ee, FindNewRowptrSize(t));
    type_ = t;
}

template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType> *SparsePartition<IndexType, ValueType>::
ExtractWindow(IndexType rs, IndexType length)
{
    if (rs + length > rowptr_size_ - 1)
        length = rowptr_size_ - rs - 1;

    SparsePartition<IndexType, ValueType> *ret = new
        SparsePartition<IndexType, ValueType>;
    IndexType es = rowptr_[rs];
    IndexType ee = rowptr_[rs+length];

    if (es == ee) return ret;
    assert(es < ee);
    iterator p_start = begin(rs);
    iterator p_end = end(rs + length - 1);
    ret->SetElems(p_start, p_end, rs + 1, ee-es, ee-es, length);
    ret->nr_rows_ = ret->rowptr_size_ - 1;
    ret->nr_cols_ = nr_cols_;
    ret->nr_nzeros_ = ret->elems_size_;
    ret->row_start_ = row_start_ + rs;
    ret->type_ = type_;

    assert(ret->nr_nzeros_ == ret->elems_.size());
    assert((size_t)(ret->rowptr_[ret->rowptr_size_-1]) == ret->elems_size_);
    return ret;
}

template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType> *SparsePartition<IndexType, ValueType>::
GetWindow(IndexType rs, IndexType length)
{
    if (rs + length > rowptr_size_ - 1)
        length = rowptr_size_ - rs - 1;

    SparsePartition<IndexType, ValueType> *ret = new
        SparsePartition<IndexType, ValueType>;
    IndexType es = rowptr_[rs];
    IndexType ee = rowptr_[rs+length];
    if (es == ee) return ret;

    // Move element range corresponding to this window
    assert(es < ee);
    ret->elems_.reserve(ee - es);
    move(elems_.begin() + es, elems_.begin() + ee, back_inserter(ret->elems_));
    ret->elems_size_ = ee - es;

    // Adjust elements' row
    typename vector<Element<IndexType, ValueType> >::iterator p_start, p_end;
    p_start = ret->elems_.begin();
    p_end = ret->elems_.end();
    for (; p_start != p_end; ++p_start) {
        (*p_start).Transform(make_pair((*p_start).GetRow() - rs,
                                       (*p_start).GetCol()));
    }

    p_start = ret->elems_.begin();
    ret->SetRowptr(p_start, p_end, length);
    ret->nr_rows_ = length;
    ret->nr_cols_ = nr_cols_;
    ret->nr_nzeros_ = ret->elems_size_;
    ret->row_start_ = row_start_ + rs;
    ret->type_ = type_;
    ret->elems_mapped_ = true;

    assert((size_t)(ret->rowptr_[ret->rowptr_size_-1]) == ret->elems_size_);
    return ret;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::
PutWindow(SparsePartition<IndexType, ValueType> *window)
{
    assert(window);
    assert(type_ == window->type_);

    IndexType rs = window->row_start_ - row_start_;
    IndexType es = rowptr_[rs];

    if (type_ == Encoding::Horizontal) {
        // Adjust element rows if putting window back to an horizontal matrix
        for (size_t i = 0; i < window->elems_size_; ++i) {
            Element<IndexType, ValueType> &e = window->elems_[i];
            e.Transform(
                make_pair(static_cast<IndexType>(e.GetRow() + rs), e.GetCol()));
        }
    }

    // Move back elements
    move(window->elems_.begin(), window->elems_.end(), elems_.begin() + es);
}

template<typename IndexType, typename ValueType>
ostream &operator<<(ostream &os, const SparsePartition<IndexType,
                    ValueType> &mat)
{
    mat.Print(os);
    return os;
}

/*
 * Builder class implemetation
 */
template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType>::
Builder::Builder(SparsePartition<IndexType, ValueType> *sp, size_t nr_rows)
    : sp_(sp),
      da_rowptr_(DynamicIndexArray(nr_rows+1))
{
    da_rowptr_.Append(0);
}

template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType>::
Builder::~Builder()
{ }

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::
Builder::NewRow(IndexType rdiff, IndexType value)
{
    for (IndexType i = 0; i < rdiff; ++i)
        da_rowptr_.Append(value);
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::
Builder::Finalize(size_t size)
{
    const IndexType &last_rowptr = da_rowptr_.GetLast();
    if (((size_t) last_rowptr) != size) {
        NewRow(1, size);
    }

    // Free old data structures
    if (sp_->rowptr_) {
        da_rowptr_.GetAllocator().destroy(sp_->rowptr_, sp_->rowptr_size_);
        da_rowptr_.GetAllocator().deallocate(sp_->rowptr_, sp_->rowptr_size_);
    }

    sp_->rowptr_size_ = da_rowptr_.GetSize();
    sp_->rowptr_ = da_rowptr_.TakeElems();
}

/* 
 * SparsePartition::iterator implemetation
 */
template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType>::iterator::iterator(
    SparsePartition *sp, IndexType r_idx)
    : sp_(sp),
      row_idx_(r_idx),
      elem_idx_(sp_->rowptr_[r_idx])
{
    assert(r_idx < sp_->rowptr_size_);
}

template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType>::iterator::iterator(
    SparsePartition *sp, IndexType r_idx, struct end_tag)
    : sp_(sp),
      row_idx_(r_idx),
      elem_idx_(sp_->rowptr_[r_idx+1])
{
    assert(r_idx < sp_->rowptr_size_);
}

template<typename IndexType, typename ValueType>
bool SparsePartition<IndexType, ValueType>::iterator::operator==(
    const iterator &pi)
{
    return !(*this != pi);
}

template<typename IndexType, typename ValueType>
bool SparsePartition<IndexType, ValueType>::iterator::operator!=(
    const iterator &pi)
{
    if (sp_ != pi.sp_)
        return true;

    if (row_idx_ != pi.row_idx_)
        return true;

    if (elem_idx_ != pi.elem_idx_)
        return true;

    return false;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::iterator::operator++()
{
    IndexType *rp = sp_->rowptr_;
    IndexType rp_size = sp_->rowptr_size_;

    assert((size_t) elem_idx_ < sp_->elems_size_);
    assert((size_t) row_idx_ < rp_size);
    while (row_idx_+1 < rp_size && rp[row_idx_+1] == elem_idx_)
        row_idx_++;

    elem_idx_++;
}

template<typename IndexType, typename ValueType>
Element<IndexType, ValueType> &
SparsePartition<IndexType, ValueType>::iterator::operator*()
{
    return sp_->elems_[elem_idx_];
}

/*
 * SparsePartitionSym class implementation
 */
template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::DivideMatrix()
{
    SparsePartition<IndexType, ValueType> *matrix = lower_matrix_;
    IndexType row_start = matrix->row_start_;
    IndexType nr_rows = matrix->rowptr_size_ - 1;
    IndexType nr_cols = matrix->nr_cols_;
    IndexType rows1 = 0;
    IndexType rows2 = 0;
    size_t elem_cnt1 = 0;
    size_t elem_cnt2 = 0;

    typename SparsePartition<IndexType, ValueType>::Builder 
        Bld1(m1_, nr_rows + 1);
    typename SparsePartition<IndexType, ValueType>::Builder
        Bld2(m2_, nr_rows + 1);
    
    m1_->type_ = Encoding::Horizontal;
    m1_->row_start_ = row_start;
    m1_->nr_cols_ = nr_cols;

    m2_->type_ = Encoding::Horizontal;
    m2_->row_start_ = row_start;
    m2_->nr_cols_ = nr_cols;
    
    for (IndexType i = 0; i < nr_rows; i++) {
        for (IndexType j = matrix->rowptr_[i]; j < matrix->rowptr_[i+1]; j++) {
            Element<IndexType, ValueType> elem(matrix->operator[](j));
            if (elem.GetCol() < row_start + 1) {
                if (rows1 < i) {
                    Bld1.NewRow(i - rows1, elem_cnt1);
                    rows1 = i;
                }

                m1_->nr_nzeros_++;
                m1_->elems_.push_back(move(elem));
                elem_cnt1++;
            } else {
                if (rows2 < i) {
                    Bld2.NewRow(i - rows2, elem_cnt2);
                    rows2 = i;
                }

                m2_->nr_nzeros_++;
                m2_->elems_.push_back(move(elem));
                elem_cnt2++;
            }
        }
    }
    
    assert(elem_cnt1 == m1_->elems_.size());
    assert(elem_cnt2 == m2_->elems_.size());
    Bld1.Finalize(elem_cnt1);
    Bld2.Finalize(elem_cnt2);
    m1_->elems_size_ = m1_->elems_.size();
    m2_->elems_size_ = m2_->elems_.size();
    m1_->nr_rows_ = m1_->rowptr_size_ - 1;
    m2_->nr_rows_ = m2_->rowptr_size_ - 1;
}

template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::MergeMatrix()
{
    SparsePartition<IndexType, ValueType> *m1 = GetFirstMatrix();
    SparsePartition<IndexType, ValueType> *m2 = GetSecondMatrix();
    SparsePartition<IndexType, ValueType> *temp = new
        SparsePartition<IndexType, ValueType>;

    SparsePartition<IndexType, ValueType> *matrix = lower_matrix_;
    IndexType row_start = matrix->row_start_;
    size_t nr_rows = matrix->rowptr_size_ - 1;
    IndexType nr_cols = matrix->nr_cols_;
    size_t elem_cnt = 0;

    typename SparsePartition<IndexType, ValueType>::Builder 
        Bld(temp, nr_rows + 1);

    temp->type_ = Encoding::Horizontal;
    temp->row_start_ = row_start;
    temp->nr_cols_ = nr_cols;
    temp->nr_nzeros_ = matrix->nr_nzeros_;
    
    for (size_t i = 0; i < nr_rows; i++) {
        if (((size_t)m1_->rowptr_size_ - 1) > i) {
            for (IndexType j = m1_->rowptr_[i]; j < m1_->rowptr_[i+1]; j++) {
                temp->elems_.push_back(m1->operator[](j));
                elem_cnt++;
            }
        }

        if (((size_t)m2_->rowptr_size_ - 1) > i) {
            for (IndexType j = m2_->rowptr_[i]; j < m2_->rowptr_[i+1]; j++) {
                temp->elems_.push_back(m2->operator[](j));
                elem_cnt++;
            }
        }
        Bld.NewRow(1, elem_cnt);
    }
    
    Bld.Finalize(elem_cnt);
    temp->nr_rows_ = temp->rowptr_size_ - 1;

    delete lower_matrix_;
    lower_matrix_ = temp;

    delete m1_;
    delete m2_;
}

template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::PrintDiagElems(ostream &out)
{
    IndexType row_start = lower_matrix_->row_start_;
    
    for (size_t i = 0; i < diagonal_size_; i++)
        cout << row_start + i + 1 << " " << row_start + i + 1 << " "
                  << diagonal_[i] << " cnt:" << i << endl;
}

template<typename IndexType, typename ValueType> template<typename IterT>
IndexType SparsePartitionSym<IndexType, ValueType>::
SetElems(IterT &pi, const IterT &pnts_end, IndexType row_start, 
         size_t limit, size_t nr_elems, size_t nr_rows)
{
    typename SparsePartition<IndexType, ValueType>::Builder 
        Bld(lower_matrix_, nr_rows);
    IndexType row_prev = 1;
    size_t elem_cnt = 0, diag_cnt = 0;

    lower_matrix_->elems_.reserve(nr_elems);
    // lower_matrix_->elems_.reserve(limit);
    for (; pi != pnts_end; ++pi) {
        Element<IndexType, ValueType> elem(*pi);
        IndexType row = elem.GetRow() - row_start + 1;
        IndexType col = elem.GetCol();
        if (row_start + row - 1 > col) {
            if (row != row_prev) {
                assert(row > row_prev);
                if (limit && diag_cnt + elem_cnt >= limit &&
                    row_prev == row - 1)
                    break;
                Bld.NewRow(row - row_prev, elem_cnt);
                row_prev = row;
            }

            // Element's row must be set to the new value
            elem.Transform(make_pair(row, col));
            lower_matrix_->elems_.push_back(move(elem));
            elem_cnt++;
        } else if (row_start + row - 1 == col) {
            diagonal_.push_back(elem.GetValue());
            diag_cnt++;
        }
    }

    // elems_'s size will never increase henceforth, so adjust capacity
    lower_matrix_->elems_.shrink_to_fit();
    lower_matrix_->elems_size_ = lower_matrix_->elems_.size();
    Bld.Finalize(lower_matrix_->elems_size_);
    diagonal_size_ = diag_cnt;
    return lower_matrix_->elems_size_ + diagonal_size_;
}

}   // end of namespace csx
}   // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_SPARSE_PARTITION_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
