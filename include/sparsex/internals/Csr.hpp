/*
 * Csr.hpp -- Routines for converting standard CSR to our intermediate
 *            representation.
 *
 * Copyright (C) 2011-2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011,      Vasileios Karakasis
 * Copyright (C) 2012-2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSR_HPP
#define SPARSEX_INTERNALS_CSR_HPP

#include "sparsex/internals/Element.hpp"
#include <inttypes.h>
#include <iterator>
#include <vector>

using namespace std;

namespace csx {

/**
 *  CSR wrapper.
 */ 
template<typename IndexType, typename ValueType>
class CSR
{
public:
    typedef IndexType idx_t;
    typedef ValueType val_t;

    CSR(const char *filename) {}  // Dummy
    CSR()
        : rowptr_(0),
          colind_(0),
          values_(0),
          nr_rows_(0),
          nr_cols_(0),
          zero_based_(false),
          symmetric_(false),
          reordered_(false)
    {}

    CSR(IndexType *rowptr, IndexType *colind, ValueType *values,
        IndexType nr_rows, IndexType nr_cols, bool zero_based)
        : rowptr_(rowptr),
          colind_(colind),
          values_(values),
          nr_rows_(nr_rows),
          nr_cols_(nr_cols),
          zero_based_(zero_based),
          symmetric_(false),
          reordered_(false)
    {
        nr_nzeros_ = rowptr_[nr_rows] - !zero_based_;
    }

    CSR(IndexType *rowptr, IndexType *colind, ValueType *values,
        IndexType nr_rows, IndexType nr_cols, bool zero_based,
        bool symmetric)
        : rowptr_(rowptr),
          colind_(colind),
          values_(values),
          nr_rows_(nr_rows),
          nr_cols_(nr_cols),
          zero_based_(zero_based),
          symmetric_(symmetric),
          reordered_(false)
    {
        nr_nzeros_ = rowptr_[nr_rows] - !zero_based_;
    }

    size_t GetNrNonzeros() const
    {
        return nr_nzeros_;
    }

    size_t GetNrRows() const 
    {
        return nr_rows_;
    }

    size_t GetNrCols() const
    {
        return nr_cols_;
    }

    bool IsSymmetric() const
    {
        return symmetric_;
    }

    bool IsZeroBased() const
    {
        return zero_based_;
    }

    bool IsReordered() const
    {
        return reordered_;
    }

    void SetReordered(vector<size_t> &inv_perm)
    {
        reordered_ = true;
        swap(permutation_, inv_perm);
    }

    // Element iterator
    class iterator;

    iterator begin()
    {
        if (reordered_)
            return iterator(this, 0, permutation_);
        else
            return iterator(this, 0);
    }

    iterator end()
    {
        if (reordered_)
            return iterator(this, nr_nzeros_, permutation_);
        else 
            return iterator(this, nr_nzeros_);
    }

    ValueType GetValue(IndexType row, IndexType col) const;
    bool SetValue(IndexType row, IndexType col, ValueType new_val);

    template<typename I, typename V>
    friend ostream &operator<<(ostream &, const CSR<I, V> &);

    IndexType *rowptr_;
    IndexType *colind_;
    ValueType *values_;
    ValueType *dvalues_;
    size_t nr_rows_, nr_cols_, nr_nzeros_;
private:
    bool zero_based_, symmetric_, reordered_;
    vector<size_t> permutation_;
};

/*
 * Returns the value of the corresponding element (0 if element doesn't exist)
 * Indexes row/col are assumed to be 1-based.
 */
// TODO: do binary search on column
template<typename IndexType, typename ValueType>
ValueType CSR<IndexType, ValueType>::GetValue(IndexType row, IndexType col) const
{
    for (IndexType i = rowptr_[row - 1] - !zero_based_;
         i < rowptr_[row] - !zero_based_; i++) {
        IndexType j = colind_[i] - !zero_based_;
        if ((j + 1) == col)
            return values_[i];
    }
    return 0;
}

/* 
 * Updates the value of the corresponding element, only if it already exists.
 * Indexes row/col are assumed to be 1-based.
 * Returns true on success, or false on error. 
 */
// TODO: do binary search on column
template<typename IndexType, typename ValueType>
bool CSR<IndexType, ValueType>::SetValue(IndexType row, IndexType col,
                                         ValueType new_val)
{
    for (IndexType i = rowptr_[row - 1] - !zero_based_;
         i < rowptr_[row] - !zero_based_; i++) {
        IndexType j = colind_[i] - !zero_based_;
        if ((j + 1) == col) {
            values_[i] = new_val;
            return true;
        }
    }
    return false;
}

template<typename IndexType, typename ValueType>
ostream &operator<<(ostream &os, const CSR<IndexType, ValueType> &mat)
{
    if (!mat.reordered_) {
        os << "Number of rows: " << mat.nr_rows_ << endl
           << "Number of cols: " << mat.nr_cols_ << endl
           << "Number of nonzeros: " << mat.nr_nzeros_ << endl;

        os << "Dump of rowptr: " << endl;
        os << "[ ";
        for (size_t i = 0; i <= mat.nr_rows_; ++i)
            os << mat.rowptr_[i] << " ";
        os << " ]" << endl;

        os << "Dump of colind: " << endl;
        os << "[ ";
        for (size_t i = 0; i < mat.nr_nzeros_; ++i)
            os << mat.colind_[i] << " ";
        os << " ]" << endl;
        
        os << "Dump of values: " << endl;
        os << "[ ";
        for (size_t i = 0; i < mat.nr_nzeros_; ++i)
            os << mat.values_[i] << " ";
        os << " ]";
    } else {
        typename CSR<IndexType, ValueType>::iterator iter = mat.begin();
        typename CSR<IndexType, ValueType>::iterator iter_end = mat.end();
            
        os << "Elements of Matrix" << endl;
        os << "------------------" << endl;
        for (; iter != iter_end; ++iter)
            os << *iter << "\n";
    }

    return os;
}

/**
 *  CSR matrix iterator.
 */ 
template<typename IndexType, typename ValueType>
class CSR<IndexType, ValueType>::iterator :
        public std::iterator<std::forward_iterator_tag,
                             Element<IndexType, ValueType> >
{
public:
    iterator(CSR<IndexType, ValueType> *csr)
        : csr_(csr),
          curr_row_(0),
          row_index_(0),
          elem_count_(0),
          curr_elem_(0)
    { }

    iterator(CSR<IndexType, ValueType> *csr, IndexType cnt)
        : csr_(csr),
          curr_row_(0),
          row_index_(0),
          elem_count_(0),
          curr_elem_(0)
    {
        for (IndexType i = 0; i < cnt; i++)
            ++(*this);
    }

    iterator(CSR<IndexType, ValueType> *csr, const vector<size_t> &inv_perm)
        : csr_(csr),
          row_index_(0),
          elem_count_(0),
          permutation_(inv_perm)
    {
        curr_row_ = permutation_[0];
        // in case first rows are empty
        while (csr_->rowptr_[curr_row_ + 1] == csr_->rowptr_[curr_row_])
            curr_row_ = permutation_[++row_index_];
        curr_elem_ = csr->rowptr_[curr_row_] - !csr->zero_based_;
    }

    iterator(CSR<IndexType, ValueType> *csr, IndexType cnt,
             const vector<size_t> &inv_perm)
        : csr_(csr),
          row_index_(0),
          elem_count_(0),
          permutation_(inv_perm)
    {
        curr_row_ = permutation_[0];
        // in case first rows are empty
        while (csr_->rowptr_[curr_row_ + 1] == csr_->rowptr_[curr_row_])
            curr_row_ = permutation_[++row_index_];
        curr_elem_ = csr->rowptr_[curr_row_] - !csr->zero_based_;
        for (IndexType i = 0; i < cnt; i++)
            ++(*this);
    }

    bool operator==(const iterator &rhs) const
    {
        return (csr_ == rhs.csr_ &&
                curr_row_ == rhs.curr_row_ &&
                curr_elem_ == rhs.curr_elem_ &&
                row_index_ == rhs.row_index_ &&
                elem_count_ == rhs.elem_count_);
    }

    bool operator!=(const iterator &rhs) const
    {
        return !(*this == rhs);
    }

    void operator++()
    {
        IndexType new_row_elem =
            csr_->rowptr_[curr_row_+1] - !csr_->zero_based_;
        ++curr_elem_;
        if  (!csr_->reordered_) {
            while (new_row_elem == curr_elem_) {
                ++curr_row_;
                new_row_elem = csr_->rowptr_[curr_row_+1] - !csr_->zero_based_;
            }

            if (curr_row_ == csr_->nr_rows_)
                assert(static_cast<size_t>(curr_elem_) == csr_->nr_nzeros_);
        } else {
            elem_count_++;
            while (new_row_elem == curr_elem_) {
                if (row_index_ < csr_->nr_rows_-1) {
                    curr_row_ = permutation_[++row_index_];
                    curr_elem_ = csr_->rowptr_[curr_row_] - !csr_->zero_based_;
                    new_row_elem = csr_->rowptr_[curr_row_+1] -
                        !csr_->zero_based_;
                }
                else {
                    curr_row_ = csr_->nr_rows_;
                    curr_elem_ = static_cast<IndexType>(csr_->nr_nzeros_);
                    break;
                }
            }

            if (row_index_ == csr_->nr_rows_)
                assert(elem_count_ == csr_->nr_nzeros_);
        }
    }

    Element<IndexType, ValueType> operator*()
    {
        assert(curr_row_ <= csr_->nr_rows_ && "out of bounds");
        assert(static_cast<size_t>(curr_elem_) < csr_->nr_nzeros_
               && "out of bounds");

        IndexType row, col;
        ValueType val;

        // Element's are one-based!
        if (csr_->reordered_) {
            row = row_index_ + 1;
        } else {
            row = curr_row_ + 1;
        }

        col = csr_->colind_[curr_elem_] + csr_->zero_based_;
        val = csr_->values_[curr_elem_];
        return Element<IndexType, ValueType>(row, col, val);
    }

private:
    CSR<IndexType, ValueType> *csr_;
    size_t curr_row_, row_index_, elem_count_;
    IndexType curr_elem_;
    vector<size_t> permutation_;
};

} // namespace csx

#endif  // SPARSEX_INTERNALS_CSR_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
