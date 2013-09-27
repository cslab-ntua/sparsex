/* -*- C++ -*-
 * 
 * csx_bench.cc -- Benchmark utilities.
 *
 * Copyright (C) 2011-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2011-2012, Vasileios Karakasis
 * Copyright (C) 2011-2012, Theodoros Gkountouvas
 * Copyright (C) 2013,      Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "csx_bench.h"
#include "csr.h"

using namespace csx;

static int GetOptionOuterLoops()
{
    const char *loops_env = getenv("OUTER_LOOPS");
    int ret = 1;
    
    if (loops_env) {
        ret = atoi(loops_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
}

static void MMFtoCSR(const char *filename, uint64_t **rowptr, uint64_t **colind,
                     double **values, size_t *nrows, size_t *ncols, size_t *nnz)
{
    MMF<uint64_t, double> mmf(filename);
    *nrows = mmf.GetNrRows(); *ncols = mmf.GetNrCols();
    *nnz = mmf.GetNrNonzeros();
	*values = (double *) malloc(sizeof(double) * mmf.GetNrNonzeros());
	*colind = (uint64_t *) malloc(sizeof(uint64_t) * mmf.GetNrNonzeros());
	*rowptr = (uint64_t *) malloc(sizeof(uint64_t) * (mmf.GetNrRows() + 1));

    MMF<uint64_t, double>::iterator iter = mmf.begin();
    MMF<uint64_t, double>::iterator iter_end = mmf.end();   
    uint64_t row_i = 0, val_i = 0, row_prev = 0;
    uint64_t row, col;
    double val;

	(*rowptr)[row_i++] = val_i;
    for (;iter != iter_end; ++iter) {
        row = (*iter).row - 1;
        col = (*iter).col - 1;
        val = (*iter).val;
		assert(row >= row_prev);
		if (row != row_prev) {
			for (uint64_t i = 0; i < row - row_prev; i++) {
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

void CheckLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    CSR<uint64_t, double> *csr = new CSR<uint64_t, double>();

    MMFtoCSR(mmf_name, &csr->rowptr_, &csr->colind_, &csr->values_, &csr->nr_rows_, &csr->nr_cols_, &csr->nr_nzeros_);
    std::cout << "Checking ... " << std::flush;
    if (!spm_mt->symmetric)
        SPMV_CHECK_FN(csr, spm_mt, 1, csr->GetNrRows(), csr->GetNrCols());
    else
        SPMV_CHECK_SYM_FN(csr, spm_mt, 1, csr->GetNrRows(), csr->GetNrCols());
    std::cout << "Check Passed" << std::endl;
    delete csr; //leak
}

void BenchLoop(spm_mt_t *spm_mt, char *mmf_name)
{
    uint64_t nrows, ncols, nnz;
    double secs, flops;
    long loops_nr = 128;

    ReadMmfSizeLine(mmf_name, nrows, ncols, nnz);
    int nr_outer_loops = GetOptionOuterLoops();
    
    for (int i = 0; i < nr_outer_loops; ++i) {
        if (!spm_mt->symmetric) {
            secs = SPMV_BENCH_FN(spm_mt, loops_nr, nrows, ncols);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            printf("m:%s f:%s s:%lu t:%lf r:%lf\n", "csx",
                   basename(mmf_name), CsxSize<double>(spm_mt), secs,
                   flops);
        } else {
            secs = SPMV_BENCH_SYM_FN(spm_mt, loops_nr, nrows, ncols);
            flops = (double)(loops_nr*nnz*2)/((double)1000*1000*secs);
            printf("m:%s f:%s ms:%lu s:%lu t:%lf r:%lf\n", "csx-sym",
                   basename(mmf_name), MapSize(spm_mt),
                   CsxSymSize<double>(spm_mt), secs, flops);
        }
    }
}
