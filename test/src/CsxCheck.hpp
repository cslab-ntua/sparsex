/*
 * Copyright (C) 2011-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/*
 * \file CsxCheck.hpp
 * \brief Checking utilities
 *
 * \author Vasileios Karakasis
 * \author Theodoros Gkountouvas
 * \author Athena Elafrou
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_TEST_CHECK_HPP
#define SPARSEX_TEST_CHECK_HPP

#include <sparsex/internals/SpmMt.hpp>

#ifdef __cplusplus
// C++ only

#include <sparsex/internals/Csr.hpp>
#include <sparsex/internals/Mmf.hpp>
#include <sparsex/internals/Types.hpp>
#include <sparsex/internals/Vector.hpp>

using namespace std;
using namespace sparsex::io;

/**
 *  Read an MMF file and create the corresponding CSR format.
 */
template<typename IndexType, typename ValueType>
void MMFtoCSR(const char *filename, IndexType **rowptr, IndexType **colind,
              ValueType **values, size_t *nrows, size_t *ncols, size_t *nnz)
{
  MMF<IndexType, ValueType> mmf(filename);
  *nrows = mmf.GetNrRows();
  *ncols = mmf.GetNrCols();
  *nnz = mmf.GetNrNonzeros();
  *values = new ValueType[mmf.GetNrNonzeros()];
  *colind = new IndexType[mmf.GetNrNonzeros()];
  *rowptr = new IndexType[mmf.GetNrRows() + 1];

  typename MMF<IndexType, ValueType>::iterator iter = mmf.begin();
  typename MMF<IndexType, ValueType>::iterator iter_end = mmf.end();   
  IndexType row_i = 0, val_i = 0, row_prev = 0;
  IndexType row, col;
  ValueType val;

  (*rowptr)[row_i++] = val_i;
  for (;iter != iter_end; ++iter) {
    row = (*iter).GetRow() - 1;
    col = (*iter).GetCol() - 1;
    val = (*iter).GetValue();
    assert(row >= row_prev);
    if (row != row_prev) {
      for (IndexType i = 0; i < row - row_prev; i++) {
	(*rowptr)[row_i++] = val_i;
      }
      row_prev = row;
    }
    (*values)[val_i] = val;
    (*colind)[val_i] = col;
    val_i++;
  }

  (*rowptr)[row_i++] = val_i;
}

/**
 *  Check the CSX SpMV result against the baseline single-thread CSR
 *  implementation.
 */
template<typename IndexType, typename ValueType>
void CheckResult(vector_t *result, ValueType alpha, vector_t *x,
                 char *matrix_file)
{
  CSR<IndexType, ValueType> *csr = new CSR<IndexType, ValueType>;
  MMFtoCSR<IndexType, ValueType>(matrix_file, &csr->rowptr_, &csr->colind_,
				 &csr->values_, &csr->nr_rows_,
				 &csr->nr_cols_, &csr->nr_nzeros_);

  cout << "Checking... " << flush;
  vector_t *y_csr = VecCreate(csr->nr_rows_);
  csr_spmv(csr, x, y_csr);
  VecScale(y_csr, y_csr, alpha);
  if (VecCompare(y_csr, result) < 0)
    exit(1);
  cout << "Check Passed" << endl;

  // Cleanup
  delete[] csr->rowptr_;
  delete[] csr->colind_;
  delete[] csr->values_;
  delete csr;
}

#endif

SPX_BEGIN_C_DECLS__
void check_result(vector_t *result, double alpha, vector_t *x,
                  char *matrix_file);
SPX_END_C_DECLS__

#endif // SPARSEX_TEST_CHECK_HPP
