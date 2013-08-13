/* -*- C++ -*-
 *
 * MatrixLoading.h --  Functions for loading sparse matrices from
 *                     different formats (currently MMF and CSR).
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef MATRIXLOADING_H__
#define MATRIXLOADING_H__

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "csr.h"
#include "mmf.h"

using namespace std;

namespace csx {

/* Forward declarations */
template<typename IndexType, typename ValueType>
class SparseInternal;

/**
 *  Loads matrix from a MMF file specifying the number of threads.
 *
 *  @param filename name of the file where the matrix is kept.
 *  @param in       buffer from which the matrix is taken.
 *  @param mmf      handler of MMF class.
 *  @param nr       number of threads to be used.
 *  @return         SparseInternal class object with the characteristics
 *                  of the matrix.
 */
template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadMMF_mt(const char *mmf_file,
                                                 const long nr_threads = 0)
{
    MMF<IndexType, ValueType> mmf(mmf_file);
    SparseInternal<IndexType, ValueType> *ret;
    ret = SparseInternal<IndexType, ValueType>::DoLoadMatrix(mmf,
                                                             nr_threads);
    return ret;
}

template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadMMF(const char *mmf_file)
{
    SparseInternal<IndexType, ValueType> *ret;

    ret = LoadMMF_mt<IndexType, ValueType>(mmf_file, 1);
    return ret;
}

/**
 *  Loads matrix from a file specifying the number of threads.
 *
 *  @param filename name of the file that matrix is kept.
 *  @param in       buffer from which the matrix is taken.
 *  @param mmf      handler of MMF class.
 *  @param nr       number of threads to be used.
 *  @return         spmsym class object with the characteristics of the 
 *                  matrix.
 */
template<typename IndexType, typename ValueType>
class SparsePartitionSym;

template<typename IndexType, typename ValueType>
SparsePartitionSym<IndexType, ValueType> *LoadMMF_Sym_mt(const char *mmf_file,
                                                         const long nr)
{
    MMF<IndexType, ValueType> mmf(mmf_file);
    SparsePartitionSym<IndexType, ValueType> *ret;
    ret = SparsePartitionSym<IndexType, ValueType>::LoadMMF_mt(mmf, nr);
    return ret;
}

// template<typename IndexType, typename ValueType>
// SparsePartitionSym<IndexType, ValueType> *LoadMMF_Sym_mt(std::ifstream &in,
//                                                          const long nr)
// {
//     MMF<IndexType, ValueType> mmf(in);
//     return SparsePartitionSym<IndexType, ValueType>::LoadMMF_mt(mmf, nr);
// }

// template<typename IndexType, typename ValueType>
// SparsePartitionSym<IndexType, ValueType> *LoadMMF_Sym_mt(const char *mmf_file,
//                                                          const long nr)
// {
//     SparsePartitionSym<IndexType, ValueType> *ret;
//     std::ifstream in;

//     in.open(mmf_file);

//     if (in.good()) {
//         ret = LoadMMF_Sym_mt<IndexType, ValueType>(in, nr);
//     } else {
//         std::cerr << "File error!" << std::endl;
//         exit(1);
//     }

//     in.close();
//     return ret;
// }

// template<typename IndexType, typename ValueType>
// SparsePartitionSym<IndexType, ValueType> *LoadMMF_mt(MMF &mmf, const long nr)
// {
//     SparsePartitionSym<IndexType, ValueType> *ret, *spm_sym;
//     long limit, cnt, row_start, nr_nzeros, n, nnz;
//     MMF::iterator iter = mmf.begin();
//     MMF::iterator iter_end = mmf.end();

//     assert(mmf.GetNrRows() == mmf.GetNrCols());
//     nr_nzeros = (mmf.GetNrNonzeros() + mmf.GetNrCols()) / 2;
//     n = mmf.GetNrCols();
//     ret = new SparsePartitionSym<IndexType, ValueType>[nr];
//     row_start = limit = cnt = 0;
//     for (long i = 0; i < nr; ++i) {
//         spm_sym = ret + i;
//         limit = (nr_nzeros - cnt) / (nr - i);
//         nnz = spm_sym->SetElems(iter, iter_end, row_start + 1, limit,
//                                 limit + 2 * mmf.GetNrRows() - 1,
//                                 mmf.GetNrRows() + 1);
//         spm_sym->lower_matrix_->SetNrNonzeros(nnz - spm_sym->GetDiagonalSize());
//         spm_sym->lower_matrix_->SetNrRows(spm_sym->lower_matrix_->GetNrRows());
//         spm_sym->lower_matrix_->SetNrCols(n);
//         spm_sym->lower_matrix_->SetRowStart(row_start);
//         row_start += spm_sym->GetDiagonalSize();
//         spm_sym->lower_matrix_->SetType(HORIZONTAL);
//         cnt += nnz;
//     }

//     assert((uint64_t) cnt == (uint64_t) nr_nzeros);
//     return ret;
// }

/**
 *  Loads matrix from CSR format specifying the number of threads.
 *
 *  @param rowptr     array "rowptr" of CSR format.
 *  @param colind     array "colind" of CSR format.
 *  @param values     array "values" of CSR format.
 *  @param nr_rows    number of rows.
 *  @param nr_cols    number of columns.
 *  @param zero_based indexing.
 *  @param nr         number of threads to be used.
 *  @return           SparseInternal class object with the characteristics
 *                    of the matrix.
 */
template<typename IndexType, typename ValueType>
SparseInternal<IndexType, ValueType> *LoadCSR_mt(IndexType *rowptr,
                                                 IndexType *colind,
                                                 ValueType *values,
                                                 IndexType nr_rows,
                                                 IndexType nr_cols,
                                                 bool zero_based,
                                                 size_t nr_threads)
{
    CSR<IndexType, ValueType> csr(rowptr, colind, values, nr_rows, nr_cols,
                                  zero_based);
    return SparseInternal<IndexType, ValueType>::DoLoadMatrix(csr, nr_threads);
}

}  //csx namespace end

#endif // MATRIXLOADING_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
