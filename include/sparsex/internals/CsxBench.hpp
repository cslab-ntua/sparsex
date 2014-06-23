/*
 * CsxBench.hpp -- Benchmark utilities.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_CSX_BENCH_HPP
#define SPARSEX_INTERNALS_CSX_BENCH_HPP

#include "sparsex/internals/Affinity.hpp"
#include "sparsex/internals/Csr.hpp"
#include "sparsex/internals/CsxMatvec.hpp"
#include "sparsex/internals/CsxUtil.hpp"
#include "sparsex/internals/Mmf.hpp"
#include "sparsex/internals/SpmMt.hpp"
#include "sparsex/internals/SpmvMethod.hpp"
#include "sparsex/internals/ThreadPool.hpp"
#include "sparsex/internals/Types.hpp"
#include "sparsex/internals/Timer.hpp"
#include "sparsex/internals/Vector.hpp"

#include <boost/thread/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/bind.hpp>

using namespace std;
using namespace csx;

extern double csx_time;

int GetOptionOuterLoops();
float spmv_bench_mt(spm_mt_t *spm_mt, size_t loops, size_t rows_nr,
                    size_t cols_nr);
void spmv_check_mt(CSR<uindex_t, spx_value_t> *spm, spm_mt_t *spm_mt,
                   size_t loops, size_t rows_nr, size_t cols_nr);
float spmv_bench_sym_mt(spm_mt_t *spm_mt, size_t loops, size_t rows_nr,
                        size_t cols_nr);
void spmv_check_sym_mt(CSR<uindex_t, spx_value_t> *spm, spm_mt_t *spm_mt,
                       size_t loops, size_t rows_nr, size_t cols_nr);

/**
 *  Check the CSX SpMV result against the baseline single-thread CSR
 *  implementation.
 *
 *  @param the (mulithreaded) CSX matrix.
 *  @param the MMF input file.
 */
template<typename IndexType, typename ValueType>
void CheckLoop(spm_mt_t *spm_mt, char *mmf_name);

/**
 *  Run CSX SpMV and record the performance information.
 */
template<typename ValueType>
void BenchLoop(spm_mt_t *spm_mt, char *mmf_name);

/**
 *  Read an MMF file and create the corresponding CSR format.
 */
template<typename IndexType, typename ValueType>
void MMFtoCSR(const char *filename, IndexType **rowptr, IndexType **colind,
              ValueType **values, size_t *nrows, size_t *ncols, size_t *nnz);

/* Implementations */
template<typename IndexType, typename ValueType>
void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    CSR<IndexType, ValueType> *csr = new CSR<IndexType, ValueType>;
    MMFtoCSR<IndexType, ValueType>(mmf_name, &csr->rowptr_, &csr->colind_,
                                   &csr->values_, &csr->nr_rows_,
                                   &csr->nr_cols_, &csr->nr_nzeros_);

    cout << "Checking... " << flush;
    if (!spm_mt->symmetric) {
        spmv_check_mt(csr, spm_mt, 1, csr->GetNrRows(), csr->GetNrCols());
    } else{
        spmv_check_sym_mt(csr, spm_mt, 1, csr->GetNrRows(), csr->GetNrCols());
    }
    cout << "Check Passed" << endl;

    // Cleanup
    delete[] csr->rowptr_;
    delete[] csr->colind_;
    delete[] csr->values_;
    delete csr;
}

template<typename ValueType>
void BenchLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    uint64_t nrows, ncols, nnz;
    double secs, flops;
    long loops_nr = 128;

    ReadMmfSizeLine(mmf_name, nrows, ncols, nnz);
    int nr_outer_loops = GetOptionOuterLoops();
    
    for (int i = 0; i < nr_outer_loops; ++i) {
        if (!spm_mt->symmetric) {
            secs = spmv_bench_mt(spm_mt, loops_nr, nrows, ncols);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            printf("m:%s f:%s s:%lu pt:%f t:%f r:%f\n", "csx",
                   basename(mmf_name), CsxSize<ValueType>(spm_mt), csx_time,
                   secs, flops);
        } else {
            secs = spmv_bench_sym_mt(spm_mt, loops_nr, nrows, ncols);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            printf("m:%s f:%s ms:%lu s:%lu pt:%f t:%f r:%f\n", "csx-sym",
                   basename(mmf_name), MapSize(spm_mt),
                   CsxSymSize<ValueType>(spm_mt), csx_time, secs, flops);
        }
    }
}

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

#endif // SPARSEX_INTERNALS_CSX_BENCH_HPP
