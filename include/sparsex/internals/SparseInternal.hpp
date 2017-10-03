/*
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file SparseInternal.hpp
 * \brief Internal representation of sparse matrices
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_INTERNALS_SPARSE_INTERNAL_HPP
#define SPARSEX_INTERNALS_SPARSE_INTERNAL_HPP

#include <sparsex/internals/Encodings.hpp>
#include <sparsex/internals/SparsePartition.hpp>

using namespace std;

namespace sparsex {
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

      template<typename I, typename V>
      friend ostream &operator<<(ostream &, const SparseInternal
				 <SparsePartition<I, V> >&);

    private:
      size_t nr_rows_, nr_cols_, nr_nzeros_, nr_partitions_; 
      PartitionType *partitions_;
    };

    template<typename PartitionType>
    template<typename InputType>
    void SparseInternal<PartitionType>::
    BuildPartitions(InputType &input, size_t nr)
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
                                         limit, limit + 2*input.GetNrRows() - 1,
                                         input.GetNrRows() + 1);
        partition->SetNrNonzeros(nnz);
        partition->SetNrRows(partition->GetRowptrSize() - 1);
        partition->SetNrCols(input.GetNrCols());
        partition->SetRowStart(row_start);
        partition->SetType(Encoding::Horizontal);
        row_start += partition->GetNrRows();
        cnt += nnz;
      }

      // assert(cnt == nr_nzeros_);
      if (cnt != nr_nzeros_) {
        LOG_ERROR << "error in input matrix (matrix has less elements than "
	  "claimed)\n";
        exit(1);
      }
    }

    template<typename IndexType, typename ValueType>
    ostream &operator<<(ostream &os, const SparseInternal
			<SparsePartition<IndexType, ValueType> > &mat)
    {
      mat.Print(os);
      return os;
    }

  }  // end of namespace csx
}  // end of namespace sparsex

#endif // SPARSEX_INTERNALS_SPARSE_INTERNAL_HPP
