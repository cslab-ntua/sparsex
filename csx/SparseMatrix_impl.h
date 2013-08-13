/* -*- C++ -*-
 *
 * SparseMatrix_impl.h --  Implementation file.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_SPARSEMATRIX_IMPL_H__
#define LIBCSX_SPARSEMATRIX_IMPL_H__

/*
 * SparseMatrix class implementation
 */
template<class MatrixType>
SparseMatrix<MatrixType>::SparseMatrix(idx_t *rowptr, idx_t *colind,
                                       val_t *values, idx_t nr_rows,
                                       idx_t nr_cols, bool zero_based)
    : MatrixType(boost::interprocess::forward<idx_t*>(rowptr),
                 boost::interprocess::forward<idx_t*>(colind), 
                 boost::interprocess::forward<val_t*>(values),
                 boost::interprocess::forward<idx_t>(nr_rows),
                 boost::interprocess::forward<idx_t>(nr_cols),
                 boost::interprocess::forward<bool>(zero_based)),
      spi_(0),
      csx_(0)
{}

// MMF-specific constructor
template<class MatrixType>
SparseMatrix<MatrixType>::SparseMatrix(const char *filename)
    : MatrixType(boost::interprocess::forward<const char*>(filename)),
      spi_(0),
      csx_(0)
{}
 
template<class MatrixType>
inline SparseMatrix<MatrixType>::~SparseMatrix()
{
    // Release resources
    if (spi_)
        delete spi_;
}

template<class MatrixType>
typename SparseMatrix<MatrixType>::val_t SparseMatrix<MatrixType>::
GetValue(idx_t row_idx, idx_t col_idx)
{
    val_t value = 0;

    if (!csx_) {
        std::cerr << "Matrix hasn't been encoded yet!" << std::endl;
    } else {
        value = GetValueCsx<idx_t, val_t>(csx_, row_idx, col_idx);
    }
    return value;
}

template<class MatrixType>
bool SparseMatrix<MatrixType>::SetValue(idx_t row_idx, idx_t col_idx,
                                        val_t value)
{
    bool set(false);

    if (!csx_) {
        std::cerr << "Matrix hasn't been encoded yet!" << std::endl;
    } else { 
        set = SetValueCsx<idx_t, val_t>(csx_, row_idx, col_idx, value);
    }
    return set;
}

template<class MatrixType>
void SparseMatrix<MatrixType>::Reorder(vector<size_t>& perm)
{
    DoReorder_RCM(*this, perm);
}

template<class MatrixType>
spm_mt *SparseMatrix<MatrixType>::CreateCsx(const RuntimeContext& rt_config,
                                            const CsxContext& csx_config,
                                            double &pre_time)
{
    spi_ = SparseInternal<idx_t, val_t>::DoLoadMatrix
        (*this, rt_config.GetNrThreads());
    csx_ = BuildCsx<idx_t, val_t>(spi_, rt_config, csx_config, pre_time);
    delete spi_;
    spi_ = 0;
    return csx_;
}

// template<class MatrixType>
// void SparseMatrix<MatrixType>::Print(std::ostream &os)
// {
//     MatrixType::Print(os);
// }

// template<class MatrixType>
// void SparseMatrix<MatrixType>::PrintEncoded(std::ostream &os) const
// {
//     if (!csx_)
//         os << "Matrix hasn't been encoded yet!" << std::endl;
//     else 
//         PrintCsx(csx_);
// }

//void Draw(const char* filename, const int width, const int height);

#endif
