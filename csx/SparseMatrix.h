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
#ifndef LIBCSX_SPARSEMATRIX_H__
#define LIBCSX_SPARSEMATRIX_H__

#include "mmf.h"
#include "csr.h"
#include "rcm.h"
//#include "draw.h"
#include "spmv.h"
#include "runtime.h"
#include "SparseInternal.h"
#include "CsxGetSet.h"
#include <boost/interprocess/detail/move.hpp>

namespace internal {

template<class MatrixType>
struct sparse_matrix_traits {
    typedef typename MatrixType::idx_t idx_t;
    typedef typename MatrixType::val_t val_t;
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

template<class MatrixType>
class SparseMatrix : public MatrixType
{
public:
    typedef typename internal::sparse_matrix_traits<MatrixType>::idx_t idx_t;
    typedef typename internal::sparse_matrix_traits<MatrixType>::val_t val_t;

    // CSR-specific constructor
    SparseMatrix(idx_t *rowptr, idx_t *colind, val_t *values,
                 idx_t nr_rows, idx_t nr_cols, bool zero_based);
    // MMF-specific constructor
    SparseMatrix(const char* filename);
    ~SparseMatrix();

    val_t GetValue(idx_t row_idx, idx_t col_idx);
    bool SetValue(idx_t row_idx, idx_t col_idx, val_t value);
    void Reorder(vector<size_t>& perm);
    spm_mt *CreateCsx(const RuntimeContext& rt_config,
                      const CsxContext& csx_config, double &time);
//    void Print(std::ostream &os); //Doesn't work for mmf
//    void PrintEncoded(std::ostream &os) const;
//    void Draw(const char* filename, const int width, const int height);

private:
    SparseInternal<idx_t, val_t> *spi_;
    spm_mt *csx_;
    AllowInstantiation<MatrixType> instantiation;
};

/**
 *  Helper template class that generates a compile time error
 *  if SparseMatrix class is instantiated with one of the template
 *  specializations that follow.
 */
// template<class MatrixType> class AllowInstantiation {};
// template<> class AllowInstantiation<MMF<double, double> >
// { AllowInstantiation() {} };
// template<> class AllowInstantiation<CSR<double, double> >
// { AllowInstantiation() {} };
// template<> class AllowInstantiation<MMF<float, float> >
// { AllowInstantiation() {} };
// template<> class AllowInstantiation<CSR<float, float> >
// { AllowInstantiation() {} };

#endif // LIBCSX_SPARSEMATRIX_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
