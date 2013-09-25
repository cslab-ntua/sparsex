/* -*- C++ -*-
 *
 * sparse_internal.h --  Internal representation of sparse matrices.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef CSX_SPARSE_INTERNAL_H__
#define CSX_SPARSE_INTERNAL_H__

#include "sparse_util.h"
#include "sparse_partition.h"
#include "runtime.h"
#include "dynarray.h"
#include "numa.h"
#include "numa_util.h"

namespace csx {

template<typename IndexType, typename ValueType>
class SparseInternal
{
public:

    SparseInternal(IndexType nr_rows, IndexType nr_cols, IndexType nr_nzeros,
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

    IndexType GetNrRows() const
    {
        return nr_rows_;
    }

    IndexType GetNrCols() const
    {
        return nr_cols_;
    }

    IndexType GetNrNonzeros() const
    {
        return nr_nzeros_;
    }

    SparsePartition<IndexType, ValueType> *GetPartition(int id)
    {
        return partitions_ + id;
    }

    template<typename MatrixType>
    void BuildPartitions(MatrixType &mat, size_t nr);

    template<typename MatrixType>
    static SparseInternal *DoLoadMatrix(MatrixType &mat, size_t nr)
    {
        SparseInternal *ret = new SparseInternal(mat.GetNrRows(),
                                                 mat.GetNrCols(),
                                                 mat.GetNrNonzeros(),
                                                 nr);
        ret->nr_rows_ = mat.GetNrRows();
        ret->nr_cols_ = mat.GetNrCols();
        ret->nr_nzeros_ = mat.GetNrNonzeros();
        ret->BuildPartitions(mat, nr);
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
    IndexType nr_rows_, nr_cols_, nr_nzeros_, nr_partitions_; 
    SparsePartition<IndexType, ValueType> *partitions_;
};

}  //csx namespace end

/*
 * SparseInternal class implementation
 */
template<typename IndexType, typename ValueType>
template<typename MatrixType>
void SparseInternal<IndexType, ValueType>::BuildPartitions(MatrixType &mat,
                                                           size_t nr)
{
    SparsePartition<IndexType, ValueType> *partition;
    unsigned long limit;
    uint64_t row_start, cnt;

    typename MatrixType::iterator iter = mat.begin();
    typename MatrixType::iterator iter_end = mat.end();

    partitions_ = new SparsePartition<IndexType, ValueType>[nr];
    row_start = limit = cnt = 0;
    for (size_t i = 0; i < nr; ++i) {
        limit = (mat.GetNrNonzeros() - cnt) / (nr - i);
        partition = partitions_ + i;
        uint64_t nnz = partition->SetElems(iter, iter_end, row_start + 1,
                                           limit, limit + mat.GetNrRows() - 1,
                                           mat.GetNrRows() + 1);
        partition->SetNrNonzeros(nnz);
        partition->SetNrRows(partition->GetRowptrSize() - 1);
        partition->SetNrCols(mat.GetNrCols());
        partition->SetRowStart(row_start);
        partition->SetType(HORIZONTAL);
        row_start += partition->GetNrRows();
        cnt += partition->GetNrNonzeros();
    }
    assert((IndexType) cnt == nr_nzeros_);
}

#endif // CSX_SPARSE_INTERNAL_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
