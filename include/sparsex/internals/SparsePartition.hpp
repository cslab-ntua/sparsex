/*
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2014, Vasileios Karakasis
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
 * \date 2011&ndash;2014
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
        type_(Encoding::None),
        elems_(NULL),
        rowptr_(NULL),
        elems_mapped_(false),
        alloc_(MemoryAllocatorImpl::GetInstance())
    {}

    ~SparsePartition()
    {
        alloc_.Destroy(rowptr_, rowptr_size_);
        if (!elems_mapped_)
            alloc_.Destroy(elems_, elems_size_);
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

    size_t GetElemsSize() const
    {
        return elems_size_;
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
    
    void SetRowPtrSize(size_t rowptr_size)
    {
        rowptr_size_ = rowptr_size;
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
     *  @param first_row  first row of the matrix.
     *  @param limit      limit for the number of elements.
     *  @param nr_elems   size for elems allocation.
     *  @param nr_rows    size for rows allocation.
     *  @return           number of elements allocated.
     */
    template <typename IterT>
    IndexType SetElems(IterT &pnts_start, const IterT &pnts_end,
                       IndexType first_row, size_t limit = 0, 
                       size_t nr_elems = 0, size_t nr_rows = 0);

    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start point iterators start.
     *  @param pnts_end   point iterators end.
     *  @param first_row  first row of the matrix.
     *  @param limit      limit for the number of elements.
     *  @param Bld        point to a builder responsible for allocations.
     *  @return           number of elements allocated.
     */
    template <typename IterT>
    IndexType SetElems(IterT &pnts_start, const IterT &pnts_end,
                       IndexType first_row, size_t limit, 
                       SparsePartition::Builder &Bld);

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
    void PutWindow(const SparsePartition *window);
 
    // SparsePartition's  iterator
    class iterator;
    iterator begin(IndexType ridx)
    {
        return iterator(this, ridx);
    }

    iterator end(IndexType ridx)
    {
        return iterator(this, ridx, true);
    }

private:
    size_t nr_rows_, nr_cols_, nr_nzeros_;
    Encoding::Type type_;
    Element<IndexType, ValueType> *elems_;
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
    typedef DynamicArray<Element<IndexType, ValueType>,
                         reallocator<Element<IndexType, ValueType>,
                                     NumaAllocator> > DynamicElemArray;
    typedef DynamicArray<IndexType,
                         reallocator<IndexType,
                                     NumaAllocator> > DynamicIndexArray;
    typedef DynamicArray<ValueType,
                     reallocator<ValueType, NumaAllocator> > DynamicValueArray;
#else
    typedef DynamicArray<Element<IndexType, ValueType> > DynamicElemArray;
    typedef DynamicArray<IndexType> DynamicIndexArray;
    typedef DynamicArray<ValueType> DynamicValueArray;
#endif

    /**
     *  Constructs a Builder object.
     *
     *  @param sp       the matrix to be built.
     *  @param nr_elems an initial guess of the nonzero elements.
     *  @param nr_rows  an initial guess of the number of rows.
     */
    Builder(SparsePartition *sp, size_t nr_rows = 0, size_t nr_elems = 0);

    ~Builder();

    /**
     *  Appends a new element in the new matrix.
     *
     *  @e the new element to be appended; <tt>e</tt> is copied during append.
     */
    void AppendElem(const Element<IndexType, ValueType> &e);

    /**
     *  Appends a new element in the new matrix.
     *
     *  @e the new element to be appended; <tt>e</tt> is moved during append.
     */
    void AppendElem(Element<IndexType, ValueType> &&e);

    /**
     *  Counts elements already allocated in the matrix.
     *
     *  @return number of elements allocated.
     */
    size_t GetElemsCnt()
    {
        return da_elems_.GetSize();
    }

    /**
     *  Allocate elements for the next row.
     *
     *  @param rdiff number of empty rows between the current and the new row.
     */
    void NewRow(IndexType rdiff = 1);

    /**
     *  Finalizes the construction of the new matrix.
     */
    void Finalize();

public:
    SparsePartition *sp_;
    DynamicElemArray da_elems_;
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
    iterator(SparsePartition *sp, IndexType row_idx, bool end = false);
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
    ValueType *diagonal_;   
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
          diagonal_(NULL),
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
        alloc_.Destroy(diagonal_, diagonal_size_);
    };
    
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
    
    ValueType *GetDiagonal()
    {
        return diagonal_;
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
        lower_matrix_->SetNrRows(nr_rows);
    }

    void SetNrCols(size_t nr_cols)
    {
        lower_matrix_->SetNrCols(nr_cols);
    }

    void SetNrNonzeros(size_t nnz)
    {
        lower_matrix_->SetNrNonzeros(nnz - GetDiagonalSize());
    }

    void SetRowStart(IndexType row_start)
    {
        lower_matrix_->SetRowStart(row_start);
    }
    
    void SetType(Encoding::Type type)
    {
        lower_matrix_->SetType(type);
    }

    size_t GetRowptrSize() const
    {
        return lower_matrix_->GetRowptrSize();
    }

    /**
     *  Divide the matrix (m) to two submatrices (m1 and m2).
     */
    void DivideMatrix();

    /**
     *  Merge the two submatrices (m1 and m2) into one matrix (m).
     */
    void MergeMatrix();
    
    class Builder;
    
    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start   point iterators start.
     *  @param pnts_end     point iterators end.
     *  @param first_row    first row of the matrix.
     *  @param limit        limit for the number of elements.
     *  @param nr_elems     size for elems allocation.
     *  @param nr_rows      size for rows allocation.
     *  @return             number of elements allocated.
     */
    template <typename IterT>
    IndexType SetElems(IterT &pnts_start, const IterT &pnts_end,
                       IndexType first_row, size_t limit = 0,
                       size_t nr_elems = 0, size_t nr_rows = 0);

    /**
     *  Function for filling the matrix using point iterators
     *
     *  @param pnts_start   point iterators start.
     *  @param pnts_end     point iterators end.
     *  @param first_row    first row of the matrix.
     *  @param limit        limit for the number of elements.
     *  @param spm_sym_bld  point to a builder responsible for allocations.
     *  @return             number of elements allocated.
     */
    template <typename IterT>
    IndexType SetElems(IterT &pnts_start, const IterT &pnts_end,
                       IndexType first_row, size_t limit,
                       SparsePartitionSym::Builder &spm_sym_bld);
                      
    void PrintDiagElems(ostream &out = cout);
};

template<typename IndexType, typename ValueType>
class SparsePartitionSym<IndexType, ValueType>::Builder
{
public:
#if SPX_USE_NUMA
    typedef DynamicArray<ValueType,
                     reallocator<ValueType, NumaAllocator> > DynamicValueArray;
#else
    typedef DynamicArray<ValueType> DynamicValueArray;
#endif

    /**
     *  Constructs a Builder object.
     *
     *  @param spm_sym  the matrix to be built.
     *  @param nr_elems an initial guess of the nonzero elements.
     *  @param nr_rows  an initial guess of the number of rows.
     */
    Builder(SparsePartitionSym *spm_sym, size_t nr_elems = 0, 
            size_t nr_rows = 0);
    virtual ~Builder();
    
    typename SparsePartition<IndexType, ValueType>::Builder *GetLowerBuilder()
    {
        return spm_bld_;
    }
    
    void AppendDiagElem(ValueType val);
    IndexType GetDiagElemsCnt();
    size_t GetElemsCnt();
    void Finalize();
    
private:
    SparsePartitionSym<IndexType, ValueType> *spm_sym_;
    typename SparsePartition<IndexType, ValueType>::Builder *spm_bld_;
    DynamicValueArray da_dvalues_;
    size_t nr_rows_;
};


/* Helper functions */

template<typename IndexType, typename ValueType>
void TestMMF(SparsePartition<IndexType, ValueType> *sp, const char *mmf_file);

/*
 *  SparsePartition class implementation
 */ 
template<typename IndexType, typename ValueType> template<typename IterT>
IndexType SparsePartition<IndexType, ValueType>::
SetElems(IterT &pi, const IterT &pnts_end, IndexType first_row, size_t limit,
         SparsePartition::Builder &Bld)
{
    // IndexType row_prev = first_row;
    IndexType row_prev = 1;
    for (; pi != pnts_end; ++pi) {
        Element<IndexType, ValueType> elem(*pi);
        IndexType row = elem.GetRow() - first_row + 1;
        if (row != row_prev) {
            assert(row > row_prev);
            if (limit && Bld.GetElemsCnt() >= limit)
                break;
            Bld.NewRow(row - row_prev);
            row_prev = row;
        }

        // New element's row must be set to the new value
        Bld.AppendElem(TransformElement(elem, make_pair(row, elem.GetCol())));
    }

    return elems_size_;
}

template<typename IndexType, typename ValueType> template<typename IterT>
IndexType SparsePartition<IndexType, ValueType>::
SetElems(IterT &pi, const IterT &pnts_end, IndexType first_row,
         size_t limit, size_t nr_elems, size_t nr_rows)
{
    SparsePartition::Builder Bld(this, nr_rows, nr_elems);
    SetElems(pi, pnts_end, first_row, limit, Bld);
    Bld.Finalize();
    return elems_size_;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::Print(std::ostream &out)
{
    iterator p_start = begin(0);
    iterator p_end = end(GetRowptrSize() - 2);
    for (auto p = p_start; p != p_end; ++p)
        out << " " << (*p);

    out << endl;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::PrintRows(ostream &out) {
    cout << "Row Ptr: ";
    for (size_t i = 0; i < rowptr_size_; i++)
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

    nr_rows_with_patterns = GetRowptrSize() - 1;
    nr_patterns = 0;
    nr_nzeros_block = 0;
    nr_transitions = 0;
    memset(nr_xform_patterns, 0, sizeof(nr_xform_patterns));
    for (size_t i = 0; i < GetRowptrSize() - 1; i++) {
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
void SparsePartition<IndexType, ValueType>::
Transform(Encoding::Type t, IndexType rs, IndexType re)
{
    if (type_ == t)
        return;

    if (!re)
        re = rowptr_size_ - 2;

    vector<Element<IndexType, ValueType> > elems;
    typename vector<Element<IndexType, ValueType> >::iterator e0, ee, es;
    typename TransformFn<IndexType>::type xform_fn =
        GetXformFn<IndexType>(type_, t);
    iterator p = begin(rs);
    iterator pe = end(re); 

    elems.reserve(elems_size_);
    // for (iterator p = p0; p != pe; ++p) {
    for (; p != pe; ++p) {
        elems.push_back(TransformElement(*p, xform_fn((*p).GetCoordinates(),
                                                      nr_rows_, nr_cols_)));
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
        re = GetRowptrSize() - 1;
        
        e0 = elems.begin();
        ee = e0;
        for (IndexType i = k; i < re; i += k) {
            es = ee;
            ee += rowptr_[i] - rowptr_[i-k];
            sort(es, ee);
        }

        es = ee;
        ee = elems.end();
        sort(es, ee);
    } else {
        e0 = elems.begin();
        ee = elems.end();
        sort(e0, ee);
    }

    if (elems_size_)
        SetElems(e0, ee, rs + 1, 0, elems_size_, FindNewRowptrSize(t));
    elems.clear();
    type_ = t;
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

    ret->elems_ = &elems_[es];
    ret->elems_size_ = ee - es;
    ret->nr_rows_ = length;
    ret->nr_cols_ = nr_cols_;
    ret->nr_nzeros_ = ret->elems_size_;
    ret->row_start_ = row_start_ + rs;

    ret->type_ = type_;
    ret->elems_mapped_ = true;
    Element<IndexType, ValueType> *ps = &elems_[es];
    ret->SetElems(ps, &elems_[ee], rs + 1, 0, ret->nr_nzeros_, length);

    assert((size_t)(ret->rowptr_[ret->rowptr_size_-1]) == ret->elems_size_);
    return ret;
}

template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType> *SparsePartition<IndexType, ValueType>::
ExtractWindow(IndexType rs, IndexType length)
{
    vector<Element<IndexType, ValueType> > elems;
    typename vector<Element<IndexType, ValueType> >::iterator
        elem_begin, elem_end;
    SparsePartition<IndexType, ValueType> *ret = new
        SparsePartition<IndexType, ValueType>();

    if (rs + length > rowptr_size_ - 1)
        length = rowptr_size_ - rs - 1;

    IndexType es = rowptr_[rs];
    IndexType ee = rowptr_[rs+length];

    assert(es <= ee);
    elems.reserve(ee - es);
    // iterator p0 = begin(rs);
    iterator p = begin(rs);
    iterator pe = end(rs + length - 1);
    // for (iterator p = p0; p != pe; ++p)
    for (; p != pe; ++p)
        elems.push_back(*p);

    elem_begin = elems.begin();
    elem_end = elems.end();
    ret->SetElems(elem_begin, elem_end, rs + 1, elems_size_);
    elems.clear();
    ret->nr_rows_ = ret->rowptr_size_ - 1;
    ret->nr_cols_ = nr_cols_;
    ret->nr_nzeros_ = ret->elems_size_;
    ret->row_start_ = row_start_ + rs;
    ret->type_ = type_;
    ret->elems_mapped_ = true;
    assert(ret->rowptr_[ret->rowptr_size_-1] == ret->elems_size_);
    return ret;
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::
PutWindow(const SparsePartition<IndexType, ValueType> *window)
{
    assert(window);
    assert(type_ == window->type_);

    IndexType rs = window->row_start_ - row_start_;
    IndexType es = rowptr_[rs];

    if (!window->elems_mapped_)
        memcpy(&elems_[es], window->elems_,
               window->elems_size_*sizeof(Element<IndexType, ValueType>));

    if (type_ == Encoding::Horizontal) {
        // Adjust element rows if putting window back to an horizontal matrix
        for (size_t i = 0; i < window->elems_size_; ++i) {
            Element<IndexType, ValueType> &e = window->elems_[i];
            window->elems_[i] = TransformElement(
                e, make_pair(static_cast<IndexType>(e.GetRow() + rs),
                             e.GetCol()));
        }
    }
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
Builder::Builder(SparsePartition<IndexType, ValueType> *sp, size_t nr_rows,
                 size_t nr_elems)
    : sp_(sp),
      da_elems_(((sp_->elems_mapped_) ?
                 DynamicElemArray(sp_->elems_, sp_->elems_size_) :
                 DynamicElemArray(nr_elems))),
      da_rowptr_(DynamicIndexArray(nr_rows+1))
{
    da_rowptr_.Append(0);
}

template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType>::
Builder::~Builder()
{ }

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::Builder::AppendElem(
    const Element<IndexType, ValueType> &e)
{
    if (sp_->elems_mapped_)
        assert(da_elems_.GetSize() < sp_->elems_size_ &&
               "out of bounds");

    da_elems_.Append(e);
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::Builder::AppendElem(
    Element<IndexType, ValueType> &&e)
{
    if (sp_->elems_mapped_)
        assert(da_elems_.GetSize() < sp_->elems_size_ &&
               "out of bounds");

    // e is an rvalue so forward it as an rvalue in order to use
    // the Append(T&&) overload and avoid the redundant copy
    da_elems_.Append(forward<Element<IndexType, ValueType> >(e));
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::
Builder::NewRow(IndexType rdiff)
{
    size_t elems_cnt = GetElemsCnt();
    for (IndexType i = 0; i < rdiff; i++)
        da_rowptr_.Append(elems_cnt);
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::
Builder::Finalize()
{
    const IndexType &last_rowptr = da_rowptr_.GetLast();
    if (((size_t) last_rowptr) != da_elems_.GetSize()) {
        NewRow();
    }

    // free old data structures
    if (!sp_->elems_mapped_ && sp_->elems_) {
        da_elems_.GetAllocator().destroy(sp_->elems_, sp_->elems_size_);
        da_elems_.GetAllocator().deallocate(sp_->elems_, sp_->elems_size_);
    }

    if (sp_->rowptr_) {
        da_rowptr_.GetAllocator().destroy(sp_->rowptr_, sp_->rowptr_size_);
        da_rowptr_.GetAllocator().deallocate(sp_->rowptr_, sp_->rowptr_size_);
    }

    if (sp_->elems_mapped_)
        assert(sp_->elems_size_ == da_elems_.GetSize());

    sp_->elems_size_ = da_elems_.GetSize();
    
    if (!sp_->elems_mapped_)
        sp_->elems_ = da_elems_.TakeElems();

    sp_->rowptr_size_ = da_rowptr_.GetSize();
    sp_->rowptr_ = da_rowptr_.TakeElems();
}

/* 
 * SparsePartition::iterator implemetation
 */
template<typename IndexType, typename ValueType>
SparsePartition<IndexType, ValueType>::iterator::iterator(
    SparsePartition *sp, IndexType r_idx, bool end)
    : sp_(sp)
{
    IndexType *rp = sp_->rowptr_;
    
    assert(r_idx < sp_->rowptr_size_);

    row_idx_ = r_idx;
    if (end) {
        elem_idx_ = rp[r_idx+1];
    } else {
        elem_idx_ = rp[r_idx];
    }
}

template<typename IndexType, typename ValueType>
bool SparsePartition<IndexType, ValueType>::iterator::operator==(
    const iterator &pi)
{
    return (sp_ == pi.sp_) &&
        (row_idx_ == pi.row_idx_) &&
        (elem_idx_ == pi.elem_idx_);
}

template<typename IndexType, typename ValueType>
bool SparsePartition<IndexType, ValueType>::iterator::operator!=(
    const iterator &pi)
{
    return !(*this == pi);
}

template<typename IndexType, typename ValueType>
void SparsePartition<IndexType, ValueType>::iterator::operator++()
{
    IndexType *rp = sp_->rowptr_;
    size_t rp_size = sp_->rowptr_size_;
    
    assert((size_t) elem_idx_ < sp_->elems_size_);
    assert((size_t) row_idx_ < rp_size);
    while (((size_t) row_idx_ + 1) < rp_size && rp[row_idx_+1] == elem_idx_)
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
    IndexType nr_nzeros = matrix->GetNrNonzeros();
    IndexType row_start = matrix->GetRowStart();
    IndexType nr_rows = matrix->GetRowptrSize() - 1;
    IndexType nr_cols = matrix->GetNrCols();
    IndexType rows1 = 0;
    IndexType rows2 = 0;

    typename SparsePartition<IndexType, ValueType>::Builder
        spmbld1(m1_, nr_rows + 1, nr_nzeros);
    typename SparsePartition<IndexType, ValueType>::Builder
        spmbld2(m2_, nr_rows + 1, nr_nzeros);
    
    m1_->SetType(Encoding::Horizontal);
    m1_->SetRowStart(row_start);
    m1_->SetNrCols(nr_cols);
    m1_->SetNrNonzeros(0);

    m2_->SetType(Encoding::Horizontal);
    m2_->SetRowStart(row_start);
    m2_->SetNrCols(nr_cols);
    m2_->SetNrNonzeros(0);
    
    for (IndexType i = 0; i < nr_rows; i++) {
        for (IndexType j = matrix->GetRow(i); j < matrix->GetRow(i+1); j++) {
            Element<IndexType, ValueType> elem = matrix->operator[](j);
            if (elem.GetCol() < row_start + 1) {
                if (rows1 < i) {
                    spmbld1.NewRow(i - rows1);
                    rows1 = i;
                }

                m1_->SetNrNonzeros(m1_->GetNrNonzeros() + 1);
                spmbld1.AppendElem(move(elem));
            } else {
                if (rows2 < i) {
                    spmbld2.NewRow(i - rows2);
                    rows2 = i;
                }

                m2_->SetNrNonzeros(m2_->GetNrNonzeros() + 1);
                spmbld2.AppendElem(move(elem));
            }
        }
    }
    
    spmbld1.Finalize();
    spmbld2.Finalize();
    m1_->SetNrRows(m1_->GetRowptrSize() - 1);
    m2_->SetNrRows(m2_->GetRowptrSize() - 1);
}

template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::MergeMatrix()
{
    SparsePartition<IndexType, ValueType> *m1 = GetFirstMatrix();
    SparsePartition<IndexType, ValueType> *m2 = GetSecondMatrix();
    SparsePartition<IndexType, ValueType> *temp = new
        SparsePartition<IndexType, ValueType>;

    SparsePartition<IndexType, ValueType> *matrix = lower_matrix_;
    IndexType row_start = matrix->GetRowStart();
    size_t nr_rows = matrix->GetRowptrSize() - 1;
    IndexType nr_cols = matrix->GetNrCols();
    IndexType nr_nzeros = matrix->GetNrNonzeros();
    
    typename SparsePartition<IndexType, ValueType>::Builder *spmbld = new
        typename SparsePartition<IndexType, ValueType>::Builder(temp,
                                                                nr_rows + 1,
                                                                nr_nzeros);

    temp->SetType(Encoding::Horizontal);
    temp->SetRowStart(row_start);
    temp->SetNrCols(nr_cols);
    temp->SetNrNonzeros(matrix->GetNrNonzeros());
    
    for (size_t i = 0; i < nr_rows; i++) {
        if ((m1_->GetRowptrSize() - 1) > i) {
            for (IndexType j = m1_->GetRow(i); j < m1_->GetRow(i+1); j++) {
                spmbld->AppendElem(m1->operator[](j));
            }
        }

        if ((m2_->GetRowptrSize() - 1) > i) {
            for (IndexType j = m2_->GetRow(i); j < m2_->GetRow(i+1); j++) {
                spmbld->AppendElem(m2->operator[](j));
            }
        }
        spmbld->NewRow();
    }
    
    spmbld->Finalize();
    temp->SetNrRows(temp->GetRowptrSize() - 1);
    
    delete spmbld;
    delete lower_matrix_;
    lower_matrix_ = temp;
    
    delete m1_;
    delete m2_;
}

template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::PrintDiagElems(ostream &out)
{
    IndexType row_start = lower_matrix_->GetRowStart();
    
    for (size_t i = 0; i < diagonal_size_; i++)
        cout << row_start + i + 1 << " " << row_start + i + 1 << " "
                  << diagonal_[i] << " cnt:" << i << endl;
}

template<typename IndexType, typename ValueType>
SparsePartitionSym<IndexType, ValueType>::Builder::
Builder(SparsePartitionSym *spm_sym, size_t nr_elems, size_t nr_rows)
    : spm_sym_(spm_sym),
      da_dvalues_(DynamicValueArray(nr_rows)),
      nr_rows_(nr_rows)
{
    SparsePartition<IndexType, ValueType> *spm = spm_sym_->GetLowerMatrix();

    spm_bld_ = new typename SparsePartition<IndexType, ValueType>::Builder
        (spm, nr_rows, nr_elems);
}

template<typename IndexType, typename ValueType>
SparsePartitionSym<IndexType, ValueType>::Builder::~Builder()
{
    delete spm_bld_;
}

template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::Builder::AppendDiagElem(
    ValueType val)
{
    assert(da_dvalues_.GetSize() < nr_rows_ && "out of bounds");
    da_dvalues_.Append(val);
}

template<typename IndexType, typename ValueType>
IndexType SparsePartitionSym<IndexType, ValueType>::Builder::GetDiagElemsCnt()
{
    return da_dvalues_.GetSize();
}

template<typename IndexType, typename ValueType>
size_t SparsePartitionSym<IndexType, ValueType>::Builder::GetElemsCnt()
{
    return GetDiagElemsCnt() + spm_bld_->GetElemsCnt();
}

template<typename IndexType, typename ValueType>
void SparsePartitionSym<IndexType, ValueType>::Builder::Finalize()
{
    spm_bld_->Finalize();
    spm_sym_->SetDiagonalSize(GetDiagElemsCnt());
    spm_sym_->SetDiagonal(da_dvalues_.TakeElems());
}

// FIXME
template<typename IndexType, typename ValueType>
template<typename IterT>
IndexType SparsePartitionSym<IndexType, ValueType>::
SetElems(IterT &pi, const IterT &pnts_end, IndexType first_row, 
         size_t limit, SparsePartitionSym::Builder &spm_sym_bld)
{
    IndexType row_prev = 1;
    for (; pi != pnts_end; ++pi) {
        IndexType row = (*pi).GetRow() - first_row + 1;
        IndexType col = (*pi).GetCol();
        if (first_row + row - 1 > col) {
            if (row != row_prev) {;
                assert(row > row_prev);
                if (limit && spm_sym_bld.GetElemsCnt() >= limit && 
                    row_prev == row - 1)
                    break;
                spm_sym_bld.GetLowerBuilder()->NewRow(row - row_prev);
                row_prev = row;
            }

            spm_sym_bld.GetLowerBuilder()->AppendElem(
                TransformElement(*pi, make_pair(row, col)));
        } else if (first_row + row - 1 == col) {
            spm_sym_bld.AppendDiagElem((*pi).GetValue());
        }
    }

    return lower_matrix_->GetElemsSize() + GetDiagonalSize();
}

template<typename IndexType, typename ValueType>
template<typename IterT>
IndexType SparsePartitionSym<IndexType, ValueType>::
SetElems(IterT &pi, const IterT &pnts_end, IndexType first_row,
         size_t limit, size_t nr_elems, size_t nr_rows)
{
    SparsePartitionSym::Builder spm_sym_bld(this, nr_elems, nr_rows);

    SetElems(pi, pnts_end, first_row, limit, spm_sym_bld);
    spm_sym_bld.Finalize();
    return lower_matrix_->GetElemsSize() + GetDiagonalSize();
}

}   // end of namespace csx
}   // end of namespace sparsex

#endif  // SPARSEX_INTERNALS_SPARSE_PARTITION_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
