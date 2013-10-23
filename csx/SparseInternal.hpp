/* -*- C++ -*-
 *
 * SparseInternal.hpp --  Internal representation of sparse matrices.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_INTERNAL_HPP
#define SPARSE_INTERNAL_HPP

#include "dynarray.h"
#include "Encodings.hpp"
#include "Runtime.hpp"
#include "SparseUtil.hpp"
#include "SparsePartition.hpp"

namespace csx {

template<typename PartitionType>
class SparseInternal
{
public:
    SparseInternal(size_t nr_rows, size_t nr_cols, size_t nr_nzeros,
                   size_t nr_partitions)
        : nr_rows_(nr_rows),
          nr_cols_(nr_cols),
          nr_nzeros_(nr_nzeros),
          nr_partitions_(nr_partitions),
          partitions_(0)
    {}

    ~SparseInternal()
    {
        delete[] partitions_;
    }

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

    PartitionType *GetPartition(int id)
    {
        return partitions_ + id;
    }

    template<typename InputType>
    void BuildPartitions(InputType &input, size_t nr);

    template<typename InputType>
    static SparseInternal *DoLoadMatrix(InputType &input, size_t nr)
    {
        SparseInternal *ret = new SparseInternal(input.GetNrRows(),
                                                 input.GetNrCols(),
                                                 input.GetNrNonzeros(),
                                                 nr);
        ret->nr_rows_ = input.GetNrRows();
        ret->nr_cols_ = input.GetNrCols();
        ret->nr_nzeros_ = input.GetNrNonzeros();
        ret->BuildPartitions(input, nr);
        return ret;
    }

    template<typename InputType>
    static SparseInternal *DoLoadMatrixSym(InputType &input, size_t nr)
    {
        assert(input.GetNrRows() == input.GetNrCols());
        SparseInternal *ret = new SparseInternal(input.GetNrRows(),
                                                 input.GetNrCols(),
                                                 input.GetNrNonzeros(),
                                                 nr);
        ret->nr_rows_ = input.GetNrRows();
        ret->nr_cols_ = input.GetNrCols();
        ret->nr_nzeros_ = (input.GetNrNonzeros() + input.GetNrCols()) / 2;
        ret->BuildPartitions(input, nr);
        return ret;
    }

    void Print(std::ostream &os) const
    {
        os << "Elements of Encoded Matrix" << endl;
        os << "--------------------------" << endl;
        for (size_t i = 0; i < nr_partitions_; i++) {
            os << "Partition " << i << ":" << endl;
            partitions_[i].Print(os);
            os << endl;
        }
    }

private:
    size_t nr_rows_, nr_cols_, nr_nzeros_, nr_partitions_; 
    PartitionType *partitions_;
};

}  //end of csx namespace

template<typename PartitionType>
template<typename InputType>
void SparseInternal<PartitionType>::BuildPartitions(InputType &input, size_t nr)
{
    PartitionType *partition;
    unsigned long limit;
    size_t cnt;
    typename InputType::idx_t row_start;
    typename InputType::iterator iter = input.begin();
    typename InputType::iterator iter_end = input.end();

    partitions_ = new PartitionType[nr];
    row_start = limit = cnt = 0;
    for (size_t i = 0; i < nr; ++i) {
        partition = partitions_ + i;
        limit = (nr_nzeros_ - cnt) / (nr - i);
        size_t nnz = partition->SetElems(iter, iter_end, row_start + 1,
                                         limit, limit + input.GetNrRows() - 1,
                                         input.GetNrRows() + 1);
        partition->SetNrNonzeros(nnz);
        partition->SetNrRows(partition->GetRowptrSize() - 1);
        partition->SetNrCols(input.GetNrCols());
        partition->SetRowStart(row_start);
        partition->SetType(Encoding::Horizontal);
        row_start += partition->GetNrRows();
        cnt += nnz;
    }

    assert(cnt == nr_nzeros_);
}

#endif // SPARSE_INTERNAL_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
