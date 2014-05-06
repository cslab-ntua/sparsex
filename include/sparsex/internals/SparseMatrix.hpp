/*
 * SparseMatrix.hpp --  Generic representation of a sparse matrix - policy-based
 *                      design.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include "sparsex/internals/Csr.hpp"
#include "sparsex/internals/CsxBuild.hpp"
#include "sparsex/internals/CsxGetSet.hpp"
#include "sparsex/internals/CsxSaveRestore.hpp"
#include "sparsex/internals/CsxUtil.hpp"
#include "sparsex/internals/Mmf.hpp"
#include "sparsex/internals/Rcm.hpp"
#include "sparsex/internals/SparseInternal.hpp"
#include "sparsex/internals/Timer.hpp"
#include "sparsex/internals/Types.hpp"

#include <boost/interprocess/detail/move.hpp>

double internal_time, csx_time, dump_time;
 
namespace internal {

template<class MatrixType>
struct matrix_traits {
    typedef typename MatrixType::idx_t idx_t;
    typedef typename MatrixType::val_t val_t;
};

/**
 *  Helper structs that generate a compile time error if SparseMatrix
 *  class is instantiated with something other than the defined template
 *  specialization, i.e., an integral and a floating point type.
 */
template<bool integral, bool floating>
void valid_instantiation();

template<>
void valid_instantiation<true, true>() {}

template<typename index_type, typename value_type>
struct allow_instantiation
{
    allow_instantiation()
    {
        valid_instantiation<std::is_integral<index_type>::value,
                            std::is_floating_point<value_type>::value>();
    }
};

template<int T>
struct Sym
{
    enum { value = T };
};

} // end namespace internal

template<class InputPolicy>
class SparseMatrix : public InputPolicy
{
public:
    typedef typename internal::matrix_traits<InputPolicy>::idx_t idx_t;
    typedef typename internal::matrix_traits<InputPolicy>::val_t val_t;

    // CSR-specific constructor
    SparseMatrix(idx_t *rowptr, idx_t *colind, val_t *values,
                 idx_t nr_rows, idx_t nr_cols, bool zero_based)
        : InputPolicy(boost::interprocess::forward<idx_t*>(rowptr),
                      boost::interprocess::forward<idx_t*>(colind), 
                      boost::interprocess::forward<val_t*>(values),
                      boost::interprocess::forward<idx_t>(nr_rows),
                      boost::interprocess::forward<idx_t>(nr_cols),
                      boost::interprocess::forward<bool>(zero_based)),
          csx_(0) {}

    // MMF-specific constructor
    SparseMatrix(const char *filename)
        : InputPolicy(boost::interprocess::forward<const char*>(filename)),
          csx_(0) {}

    ~SparseMatrix()
    {
        // if (csx_)
        //     Destroy();
    }

    void Reorder(vector<size_t>& perm)
    {
        DoReorder_RCM(*this, perm);
    }

    spm_mt_t *CreateCsx()
    {
        spm_mt_t *spm = 0;
        RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();

        try {
            if (config.GetProperty<bool>
                (RuntimeConfiguration::MatrixSymmetric)) {
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

    val_t GetValue(idx_t row_idx, idx_t col_idx)
    {
        val_t value = 0;

        if (!csx_) {
            LOG_ERROR << "matrix hasn't been encoded yet\n";
            exit(1);
        } else {
            if (csx_->symmetric) {
                if (!GetValueCsxSym<idx_t, val_t>(csx_, row_idx,
                                                  col_idx, &value))
                    LOG_WARNING << "matrix entry not found\n";
            } else {
                if (!GetValueCsx<idx_t, val_t>(csx_, row_idx,
                                               col_idx, &value))
                    LOG_WARNING << "matrix entry not found\n";
            }
        }

        return value;
    }

    bool SetValue(idx_t row_idx, idx_t col_idx, val_t value)
    {
        if (!csx_) {
            LOG_ERROR << "matrix hasn't been encoded yet\n";
            exit(1);
        } else { 
            if (csx_->symmetric) {
                if (!SetValueCsxSym<idx_t, val_t>(csx_, row_idx,
                                                  col_idx, value)) {
                    LOG_WARNING << "matrix entry not found\n";
                    return false;
                }
            } else {
                if (!SetValueCsx<idx_t, val_t>(csx_, row_idx,
                                               col_idx, value)) {
                    LOG_WARNING << "matrix entry not found\n";
                    return false;
                }
            }
        }

        return true;
    }

    void Save(const char *filename)
    {
        timing::Timer timer;
        if (!csx_) {
            LOG_ERROR << "Matrix hasn't been encoded yet\n";
            exit(1);
        } else {
            timer.Start();
            SaveCsx<idx_t, val_t>(csx_, filename, NULL);
            timer.Pause();
            dump_time = timer.ElapsedTime();
        }
    }

    void Destroy()
    {
        PutSpmMt<val_t>(csx_);
        csx_ = 0;
    }

private:
    spm_mt_t *csx_;
    internal::allow_instantiation<idx_t, val_t> instantiation;

private:
    spm_mt_t *CreateCsx(internal::Sym<true>)
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

    spm_mt_t *CreateCsx(internal::Sym<false>)
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
};

#endif // SPARSE_MATRIX_HPP

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
