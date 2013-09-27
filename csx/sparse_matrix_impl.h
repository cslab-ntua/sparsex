/* -*- C++ -*-
 *
 * sparse_matrix_impl.h --  Implementation file.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_MATRIX_IMPL_H__
#define SPARSE_MATRIX_IMPL_H__

using namespace std;

/*
 * SparseMatrix class implementation
 */
template<class InputPolicy>
inline 
SparseMatrix<InputPolicy>::SparseMatrix(idx_t *rowptr, idx_t *colind,
                                        val_t *values, idx_t nr_rows,
                                        idx_t nr_cols, bool zero_based)
    : InputPolicy(boost::interprocess::forward<idx_t*>(rowptr),
                  boost::interprocess::forward<idx_t*>(colind), 
                  boost::interprocess::forward<val_t*>(values),
                  boost::interprocess::forward<idx_t>(nr_rows),
                  boost::interprocess::forward<idx_t>(nr_cols),
                  boost::interprocess::forward<bool>(zero_based)),
      csx_(0)
{}

template<class InputPolicy>
inline 
SparseMatrix<InputPolicy>::SparseMatrix(const char *filename)
    : InputPolicy(boost::interprocess::forward<const char*>(filename)),
      csx_(0)
{}
 
template<class InputPolicy>
inline
SparseMatrix<InputPolicy>::~SparseMatrix()
{
    // if (csx_)
    //     Destroy();
}

template<class InputPolicy>
inline 
void SparseMatrix<InputPolicy>::Reorder(vector<size_t>& perm)
{
    DoReorder_RCM(*this, perm);
}

template<class InputPolicy>
template<bool IsSym>
inline
spm_mt_t *SparseMatrix<InputPolicy>::CreateCsx(double &pre_time)
{
    spm_mt_t *spm = 0;

    try {
        spm = CreateCsx(pre_time, internal::Sym<IsSym>());
    } catch (std::exception &e) {
        LOG_ERROR << e.what() << "\n";
        exit(1);
    }

    return spm;
}

template<class InputPolicy>
spm_mt_t *SparseMatrix<InputPolicy>::CreateCsx(double &pre_time,
                                               internal::Sym<false>)
{
    RuntimeContext& rt_config = RuntimeContext::GetInstance();
    SparseInternal<idx_t, val_t> *spi;

    spi = SparseInternal<idx_t, val_t>::
        DoLoadMatrix(*this, rt_config.GetNrThreads());
    csx_ = BuildCsx<idx_t, val_t>(spi, pre_time);
    delete spi;
    spi = 0;

    return csx_;
}

template<class InputPolicy>
spm_mt_t *SparseMatrix<InputPolicy>::CreateCsx(double &pre_time, 
                                               internal::Sym<true>)
{
    RuntimeContext& rt_config = RuntimeContext::GetInstance();
    SparsePartitionSym<idx_t, val_t> *spms_sym;

    spms_sym = SparsePartitionSym<idx_t, val_t>::
        DoLoadMatrix(*this, rt_config.GetNrThreads());
    csx_ = BuildCsxSym<idx_t, val_t>(spms_sym, pre_time);
    //FIXME memory leak:this causes seg (problem somewhere in PreprocessThreadSym)
    delete[] spms_sym;
    spms_sym = 0;

    return csx_;
}

template<class InputPolicy>
inline
typename SparseMatrix<InputPolicy>::val_t SparseMatrix<InputPolicy>::
GetValue(idx_t row_idx, idx_t col_idx)
{
    val_t value = 0;

    if (!csx_) {
        LOG_ERROR << "matrix hasn't been encoded yet\n";
        exit(1);
    } else {
        if (csx_->symmetric) {
            if (!GetValueCsxSym<idx_t, val_t>(csx_, row_idx, col_idx, &value))
                LOG_WARNING << "matrix entry not found\n";
        } else {
            if (!GetValueCsx<idx_t, val_t>(csx_, row_idx, col_idx, &value))
                LOG_WARNING << "matrix entry not found\n";
        }
    }

    return value;
}

template<class InputPolicy>
inline
bool SparseMatrix<InputPolicy>::SetValue(idx_t row_idx, idx_t col_idx,
                                         val_t value)
{
    if (!csx_) {
        LOG_ERROR << "matrix hasn't been encoded yet\n";
        exit(1);
    } else { 
        if (csx_->symmetric) {
            if (!SetValueCsxSym<idx_t, val_t>(csx_, row_idx, col_idx, value)) {
                LOG_WARNING << "matrix entry not found\n";
                return false;
            }
        } else {
            if (!SetValueCsx<idx_t, val_t>(csx_, row_idx, col_idx, value)) {
                LOG_WARNING << "matrix entry not found\n";
                return false;
            }
        }
    }

    return true;
}

template<class InputPolicy>
inline
void SparseMatrix<InputPolicy>::Save(const char *filename)
{
    if (!csx_) {
        LOG_ERROR << "Matrix hasn't been encoded yet\n";
        exit(1);
    } else { 
        SaveCsx<idx_t, val_t>(csx_, filename, NULL);//FIXME
    }
}

template<class InputPolicy>
inline
void SparseMatrix<InputPolicy>::Destroy()
{
    PutSpmMt<val_t>(csx_);
    csx_ = 0;
}

#endif // SPARSE_MATRIX_IMPL_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
