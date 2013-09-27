/* -*- C++ -*-
 *
 * sparse_matrix.h --  Generic representation of a sparse matrix.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_MATRIX_H__
#define SPARSE_MATRIX_H__

#include "mmf.h"
#include "csr.h"
#include "rcm.h"
#include "csx_build.h"
#include "sparse_internal.h"
#include "csx_util.h"
#include "csx_get_set.h"
#include "csx_save_restore.h"

#include <boost/interprocess/detail/move.hpp>

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

} // end namespace internal

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
    template<bool IsSym>
    spm_mt_t *CreateCsx(double &time);
    val_t GetValue(idx_t row_idx, idx_t col_idx);
    bool SetValue(idx_t row_idx, idx_t col_idx, val_t value);
    void Save(const char *filename);
    void Destroy();

private:
    spm_mt_t *csx_;
    AllowInstantiation<InputPolicy> instantiation;

    spm_mt_t *CreateCsx(double &pre_time, internal::Sym<true>);
    spm_mt_t *CreateCsx(double &pre_time, internal::Sym<false>);
};

#endif // SPARSE_MATRIX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
