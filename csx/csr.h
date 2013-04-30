/* -*- C++ -*-
 *
 * csr.h -- Routines for converting standard CSR to our intermediate
 *          representation.
 *
 * Copyright (C) 2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSR_H__
#define CSR_H__

#include <inttypes.h>
#include <iterator>
#include "spm_bits.h"

namespace csx {

/**
 *  Read-only wrapper of CSR.
 */ 
template<typename IndexType, typename ValueType>
class CSR
{
public:
    CSR(IndexType *rowptr, IndexType *colind, ValueType *values,
        IndexType nr_rows, IndexType nr_cols, bool zero_based)
        : rowptr_(rowptr),
          colind_(colind),
          values_(values),
          nr_rows_(nr_rows),
          nr_cols_(nr_cols),
          zero_based_(zero_based),
          symmetric_(false) //FIXME: give symmetric argument
    {
        nr_nzeros_ = rowptr_[nr_rows] - !zero_based_;
    }

    uint64_t GetNrNonzeros() const
    {
        return static_cast<uint64_t>(nr_nzeros_);
    }

    uint64_t GetNrRows() const 
    {
        return static_cast<uint64_t>(nr_rows_);
    }

    uint64_t GetNrCols() const
    {
        return static_cast<uint64_t>(nr_cols_);
    }

    bool IsZeroBased() { return zero_based_; }
    bool IsSymmetric() { return symmetric_; }

    // CooElem iterator
    class iterator;

    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
        return iterator(this, nr_nzeros_);
    }

    void Print(std::ostream &os) const
    {
        os << "Number of rows: " << nr_rows_ << std::endl
           << "Number of cols: " << nr_cols_ << std::endl
           << "Number of nonzeros: " << nr_nzeros_ << std::endl;

        os << "Dump of rowptr: " << std::endl;
        os << "[ ";
        for (IndexType i = 0; i <= nr_rows_; ++i)
            os << rowptr_[i] << " ";
        os << " ]" << std::endl;

        os << "Dump of colind: " << std::endl;
        os << "[ ";
        for (IndexType i = 0; i < nr_nzeros_; ++i)
            os << colind_[i] << " ";
        os << " ]" << std::endl;
        
        os << "Dump of values: " << std::endl;
        os << "[ ";
        for (IndexType i = 0; i < nr_nzeros_; ++i)
            os << values_[i] << " ";
        os << " ]";
    }

    //ValueType GetValue(IndexType row, IndexType col);
    //bool SetValue(IndexType row, IndexType col, ValueType new_val);

//private:
    IndexType *rowptr_;
    IndexType *colind_;
    ValueType *values_;
    IndexType nr_rows_;
    IndexType nr_cols_;
    IndexType nr_nzeros_;
    bool zero_based_, symmetric_;
};

template<typename IndexType, typename ValueType>
std::ostream &operator<<(std::ostream &os,
                         const CSR<IndexType, ValueType> &mat)
{
    mat.Print(os);
    return os;
}

/**
 *  CSR matrix iterator.
 */ 
template<typename IndexType, typename ValueType>
class CSR<IndexType, ValueType>::iterator :
        public std::iterator<std::forward_iterator_tag, CooElem>
{
public:
    iterator()
    {
        assert(false && "no default iterator");
    }

    iterator(CSR<IndexType, ValueType> *csr)
        : csr_(csr),
          curr_row_(0),
          curr_elem_(0)
    { }

    iterator(CSR<IndexType, ValueType> *csr, IndexType cnt)
        : csr_(csr),
          curr_row_(0),
          curr_elem_(0)
    {
        for (IndexType i = 0; i < cnt; i++)
            ++(*this);
    }

    bool operator==(const iterator &rhs) const
    {
        return (csr_ == rhs.csr_ &&
                curr_row_ == rhs.curr_row_ &&
                curr_elem_ == rhs.curr_elem_);
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
        while (new_row_elem == curr_elem_) {
            ++curr_row_;
            new_row_elem = csr_->rowptr_[curr_row_+1] - !csr_->zero_based_;
        }

        if (curr_row_ == csr_->nr_rows_)
            assert(curr_elem_ == csr_->nr_nzeros_);
    }

    CooElem operator*()
    {
        assert(curr_row_ <= csr_->nr_rows_ && "out of bounds");
        assert(curr_elem_ < csr_->nr_nzeros_ && "out of bounds");

        CooElem ret;    // CooElem's are one-based!
        ret.y = static_cast<uint64_t>(curr_row_ + 1);
        ret.x = static_cast<uint64_t>(csr_->colind_[curr_elem_] + csr_->zero_based_);
        ret.val = static_cast<double>(csr_->values_[curr_elem_]);
//         std::cout << "(" << ret.y << ", " << ret.x << ", " << ret.val
//                   << ")" << std::endl;
        return ret;
    }

private:
    CSR<IndexType, ValueType> *csr_;
    IndexType curr_row_, curr_elem_;
};

} // namespace csx

#endif  // CSR_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
