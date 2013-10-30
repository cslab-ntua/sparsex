/* -*- C++ -*-
 *
 * SparseMatrix.hpp --  Generic representation of a sparse matrix.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include "../C-API/mattype.h"
#include "Csr.hpp"
#include "CsxBuild.hpp"
#include "CsxGetSet.hpp"
#include "CsxSaveRestore.hpp"
#include "CsxUtil.hpp"
#include "Mmf.hpp"
#include "Rcm.hpp"
#include "SparseInternal.hpp"
#include "Timer.hpp"

#include <boost/interprocess/detail/move.hpp>

double internal_time, csx_time, dump_time;
 
namespace internal {

template<class MatrixType>
struct sparse_matrix_traits {
    typedef typename MatrixType::idx_t idx_t;
    typedef typename MatrixType::val_t val_t;
};

template<int T>
struct Sym
{
    enum { value = T };
};

/**
 *  Helper template class that generates a compile time error
 *  if SparseMatrix class is instantiated with something other
 *  than this class' template specializations.
 */
template<class MatrixType> class AllowInstantiation { AllowInstantiation() {} };
template<> class AllowInstantiation<MMF<uint64_t, double> > {};
template<> class AllowInstantiation<CSR<uint64_t, double> > {};
template<> class AllowInstantiation<MMF<uint32_t, double> > {};
template<> class AllowInstantiation<CSR<uint32_t, double> > {};
template<> class AllowInstantiation<MMF<int, double> > {};
template<> class AllowInstantiation<CSR<int, double> > {};

} // end namespace internal

template<class InputPolicy>
class SparseMatrix : public InputPolicy
{
public:
    typedef typename internal::sparse_matrix_traits<InputPolicy>::idx_t idx_t;
    typedef typename internal::sparse_matrix_traits<InputPolicy>::val_t val_t;

    // CSR-specific constructor
    SparseMatrix(idx_t *rowptr, idx_t *colind, val_t *values,
                 idx_t nr_rows, idx_t nr_cols, bool zero_based);
    // MMF-specific constructor
    SparseMatrix(const char* filename);
    ~SparseMatrix();

    void Reorder(vector<size_t>& perm);
    spm_mt_t *CreateCsx();
    val_t GetValue(idx_t row_idx, idx_t col_idx);
    bool SetValue(idx_t row_idx, idx_t col_idx, val_t value);
    void Save(const char *filename);
    void Destroy();

private:
    spm_mt_t *csx_;
    internal::AllowInstantiation<InputPolicy> instantiation;

    spm_mt_t *CreateCsx(internal::Sym<true>);
    spm_mt_t *CreateCsx(internal::Sym<false>);
};

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
inline
spm_mt_t *SparseMatrix<InputPolicy>::CreateCsx()
{
    spm_mt_t *spm = 0;
    RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();

    try {
        if (config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric)) {
            spm = CreateCsx(internal::Sym<true>());
        } else {
            spm = CreateCsx(internal::Sym<false>());
        }
    } catch (std::exception &e) {
        LOG_ERROR << e.what() << "\n";
        exit(1);
    }

    return spm;
}

template<class InputPolicy>
spm_mt_t *SparseMatrix<InputPolicy>::CreateCsx(internal::Sym<false>)
{
    RuntimeContext& rt_context = RuntimeContext::GetInstance();
    SparseInternal<SparsePartition<idx_t, val_t> > *spi;
    timing::Timer timer;

    // Converting to internal representation
    timer.Start();
    spi = SparseInternal<SparsePartition<idx_t, val_t> >::
        DoLoadMatrix(*this, rt_context.GetNrThreads());
    timer.Pause();
    internal_time = timer.ElapsedTime();

    // Converting to CSX
    timer.Clear();
    timer.Start();
    csx_ = BuildCsx<idx_t, val_t>(spi);
    timer.Pause();
    csx_time = timer.ElapsedTime();

    // Clenaup
    delete spi;
    return csx_;
}

template<class InputPolicy>
spm_mt_t *SparseMatrix<InputPolicy>::CreateCsx(internal::Sym<true>)
{
    RuntimeContext& rt_context = RuntimeContext::GetInstance();
    SparseInternal<SparsePartitionSym<idx_t, val_t> > *spi;
    timing::Timer timer;
                   
    // Converting to internal representation
    timer.Start();
    spi = SparseInternal<SparsePartitionSym<idx_t, val_t> >::
        DoLoadMatrixSym(*this, rt_context.GetNrThreads());
    timer.Pause();
    internal_time = timer.ElapsedTime();

    // Converting to CSX
    timer.Clear();
    timer.Start();
    csx_ = BuildCsxSym<idx_t, val_t>(spi);
    timer.Pause();
    csx_time = timer.ElapsedTime();

    // Clenaup
    delete spi;
    spi = 0;

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
    timing::Timer timer;
    if (!csx_) {
        LOG_ERROR << "Matrix hasn't been encoded yet\n";
        exit(1);
    } else {
        timer.Start();
        SaveCsx<idx_t, val_t>(csx_, filename, NULL);//FIXME
        timer.Pause();
        dump_time = timer.ElapsedTime();
    }
}

template<class InputPolicy>
inline
void SparseMatrix<InputPolicy>::Destroy()
{
    PutSpmMt<val_t>(csx_);
    csx_ = 0;
}

extern template class SparseMatrix<MMF<index_t, value_t> >;
extern template class SparseMatrix<CSR<index_t, value_t> >;

#endif // SPARSE_MATRIX_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
