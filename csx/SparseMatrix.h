/* -*- C++ -*-
 *
 * SparseMatrix.h --  Generic representation of a sparse matrix.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSEMATRIX_H__
#define SPARSEMATRIX_H__

#include "mmf.h"
#include "csr.h"
#include "rcm.h"
#include "spmv.h"
#include "runtime.h"
#include "SparseInternal.h"
#include <boost/interprocess/detail/move.hpp>

template<class MatrixType>
struct sparse_matrix_traits {
    typedef typename MatrixType::index_t index_t;
    typedef typename MatrixType::value_t value_t;
};

template<class MatrixType>
class SparseMatrix : public MatrixType
{
public:
    typedef typename sparse_matrix_traits<MatrixType>::index_t index_t;
    typedef typename sparse_matrix_traits<MatrixType>::value_t value_t;

    // CSR-specific constructor
    SparseMatrix(index_t *rowptr, index_t *colind, value_t *values,
                 index_t nr_rows, index_t nr_cols, bool zero_based,
                 size_t nr_threads);
    // MMF-specific constructor
    SparseMatrix(const char* filename);
    ~SparseMatrix();

    void Print(std::ostream &os); //Doesn't work for mmf
    void PrintEncoded(std::ostream &os) const;
    void Reorder();
    value_t GetValue(index_t row_idx, index_t col_idx);
    void SetValue(index_t row_idx, index_t col_idx, value_t value);

    spm_mt *CreateCsx(const RuntimeContext& rt_context);

private:
    SparseInternal<index_t, value_t> *spi_;
    spm_mt *csx_;
};


/*
 * SparseMatrix class implementation
 */
template<class MatrixType>
SparseMatrix<MatrixType>::SparseMatrix(index_t *rowptr, index_t *colind,
                                       value_t *values, index_t nr_rows,
                                       index_t nr_cols, bool zero_based,
                                       size_t nr_threads)
    : MatrixType(boost::interprocess::forward<index_t*>(rowptr),
                 boost::interprocess::forward<index_t*>(colind), 
                 boost::interprocess::forward<value_t*>(values),
                 boost::interprocess::forward<index_t>(nr_rows),
                 boost::interprocess::forward<index_t>(nr_cols),
                 boost::interprocess::forward<bool>(zero_based),
                 boost::interprocess::forward<size_t>(nr_threads)),
      spi_(NULL),
      csx_(NULL)
{}

// MMF-specific constructor
template<class MatrixType>
SparseMatrix<MatrixType>::SparseMatrix(const char* filename)
    : MatrixType(boost::interprocess::forward
                 <ifstream&>(*(new ifstream(filename)))),
      spi_(NULL),
      csx_(NULL)
{}
 
template<class MatrixType>
SparseMatrix<MatrixType>::~SparseMatrix()
{
    // Release resources
}

template<class MatrixType>
void SparseMatrix<MatrixType>::Print(std::ostream &os)
{
    MatrixType::Print(os);
}

template<class MatrixType>
void SparseMatrix<MatrixType>::PrintEncoded(std::ostream &os) const
{
    if (!spi_)
        std::cout << "Matrix hasn't been encoded yet!" << std::endl;
    else 
        spi_->Print(os);
}

template<class MatrixType>
void SparseMatrix<MatrixType>::Reorder()
{
    DoReorder_RCM(*this);
}

template<class MatrixType>
spm_mt *SparseMatrix<MatrixType>::CreateCsx(const RuntimeContext& rt_context)
{
    spi_ = SparseInternal<index_t, value_t>::DoLoadMatrix
        (*this, rt_context.GetNrThreads());
    csx_ = BuildCsx<index_t, value_t>(spi_);
    return csx_;
}

template<class MatrixType>
typename SparseMatrix<MatrixType>::value_t SparseMatrix<MatrixType>::
GetValue(index_t row_idx, index_t col_idx)
{
    value_t value;

    if (!csx_) {
        // GetValue from CSR/MMF(attention when iterator is on the file)
        value = GetValue(row_idx, col_idx);
    } else {
        //value = GetValueCsx<index_t, value_t>(csx_, row_idx, col_idx);
    }
    return value;
}

template<class MatrixType>
void SparseMatrix<MatrixType>::SetValue(index_t row_idx, index_t col_idx,
                                        value_t value)
{
    //bool set;

    if (!csx_) {
        // SetValue to CSR/MMF
    } else { 
        //SetValueCsx<index_t, value_t>(csx_, row_idx, col_idx, value);
    }
}

#endif // SPARSEMATRIX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
