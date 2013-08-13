/* -*- C++ -*-
 *
 * SparseMatrixWrapper.cc -- Explicit instantiation definitions of the SparseMatrix
 *                           class and wrappers of the SparseMatrix routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "SparseMatrixWrapper.h"
#include "SparseMatrix.h"
#include "SparseMatrix_impl.h"
#include "CsxSaveRestore.h"
#include "CsxGetSet.h"
#include "CsxUtil.h"
#include "csx.h"
#include "timer.h"

/**
 *  Explicit instantiation definitions.
 */
template class SparseMatrix<CSR<index_t, value_t> >;
template class SparseMatrix<MMF<index_t, value_t> >;

template value_t GetValueCsx<index_t, value_t>(void *, index_t, index_t);
template bool SetValueCsx<index_t, value_t>(void *, index_t, index_t, value_t);
template void SaveCsx<index_t, value_t>(void *, const char *, bool, index_t *);
template spm_mt_t *RestoreCsx<index_t, value_t>(const char *,
                                                const RuntimeContext &,
                                                index_t **);
template void PutSpmMt<index_t, value_t>(spm_mt_t *);

/**
 *  Wrapper routines.
 */
extern "C" {

void *CreateCSR(index_t *rowptr, index_t *colind, value_t *values,
                index_t nr_rows, index_t nr_cols, int zero_based)
{
    SparseMatrix<CSR<index_t, value_t> > *matrix =
        new SparseMatrix<CSR<index_t, value_t> >(rowptr, colind, values,
                                                 nr_rows, nr_cols, zero_based);
    if (!matrix)
        return 0;
    return (void *)matrix;
}

void *CreateMMF(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz)
{
    SparseMatrix<MMF<index_t, value_t> > *matrix =
        new SparseMatrix<MMF<index_t, value_t> >(filename);
    *nr_rows = (uint64_t)matrix->GetNrRows();
    *nr_cols = (uint64_t)matrix->GetNrCols();
    *nnz = (uint64_t)matrix->GetNrNonzeros();

    if (!matrix)
        return NULL;
    return (void *)matrix;
}

void *ReorderCSR(void *matrix, index_t **permutation)
{
    SparseMatrix<CSR<index_t, value_t> > *mat = 
        (SparseMatrix<CSR<index_t, value_t> > *)matrix;
    *permutation = (index_t *) xmalloc(mat->GetNrRows()*sizeof(index_t));
    std::vector<size_t> perm;
    mat->Reorder(perm);
    std::copy(perm.begin(), perm.end(), *permutation);  // FIXME avoid copy
    return (void *)mat;
}

void *ReorderMMF(void *matrix, index_t **permutation) 
{
    SparseMatrix<MMF<index_t, value_t> > *mat = 
        (SparseMatrix<MMF<index_t, value_t> > *)matrix;

    std::vector<size_t> perm;
    mat->Reorder(perm);
    *permutation = (index_t *) xmalloc(mat->GetNrRows()*sizeof(index_t));
    std::copy(perm.begin(), perm.end(), *permutation);  // FIXME avoid copy

    return (void *)mat;
}

void *TuneCSR(void *matrix, double *pre_time)
{
    SparseMatrix<CSR<index_t, value_t> > *mat = 
        (SparseMatrix<CSR<index_t, value_t> > *)matrix;

    /* Load runtime configuration */
    RuntimeContext &rt_config = RuntimeContext::GetInstance();
    CsxContext csx_config;
    Configuration config;
    config = ConfigFromEnv(config, false, true);    // FIXME
    rt_config.SetRuntimeContext(config);
    csx_config.SetCsxContext(config);

    spm_mt_t *spm_mt = mat->CreateCsx(rt_config, csx_config, *pre_time);
    return (void *)spm_mt;
}

void *TuneMMF(void *matrix, double *pre_time)
{
    SparseMatrix<MMF<index_t, value_t> > *mat = 
        (SparseMatrix<MMF<index_t, value_t> > *)matrix;

    /* Load runtime configuration */
    RuntimeContext &rt_config = RuntimeContext::GetInstance();
    CsxContext csx_config;
    Configuration config;
    config = ConfigFromEnv(config, false, true);    // FIXME
    rt_config.SetRuntimeContext(config);
    csx_config.SetCsxContext(config);

    spm_mt_t *spm_mt = mat->CreateCsx(rt_config, csx_config, *pre_time);
    return (void *)spm_mt;
}

void DestroyCSR(void *matrix)
{
    SparseMatrix<CSR<index_t, value_t> > *mat = 
        (SparseMatrix<CSR<index_t, value_t> > *)matrix;
    delete mat;
}

void DestroyMMF(void *matrix)
{
    SparseMatrix<MMF<index_t, value_t> > *mat = 
        (SparseMatrix<MMF<index_t, value_t> > *)matrix;
    delete mat;
}

void SaveTuned(void *matrix, const char *filename, index_t *permutation)
{
    spm_mt_t *spm_mt = (spm_mt_t *)matrix;
    SaveCsx<index_t, value_t>(spm_mt, filename, false, permutation); //FIXME
}

void *LoadTuned(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz, int *symmetric, index_t **permutation,
                double *pre_time)
{
    RuntimeContext &rt_config = RuntimeContext::GetInstance();
    spm_mt_t *spm_mt = 0;
    csx_t<index_t, value_t> *csx = 0;
    csx::Timer timer;
    
    timer.Start();
    spm_mt = RestoreCsx<index_t, value_t>(filename, rt_config, permutation);
    timer.Pause();
    *pre_time = timer.ElapsedTime();
	for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        csx = (csx_t<index_t, value_t> *) spm_mt->spm_threads[i].spm;//cast to csx_double
        *nr_rows += csx->nrows;
        *nnz += csx->nnz;
    }
    *nr_cols = csx->ncols;
    *symmetric = spm_mt->symmetric;

    return (void *)spm_mt;
}

value_t GetValue(void *matrix, index_t row, index_t col)
{
    spm_mt_t *spm_mt = (spm_mt_t *)matrix;
    value_t value = GetValueCsx<index_t, value_t>(spm_mt, row, col);
    return value;
}

int SetValue(void *matrix, index_t row, index_t col, value_t value)
{
    spm_mt_t *spm_mt = (spm_mt_t *)matrix;
    if (SetValueCsx<index_t, value_t>(spm_mt, row, col, value))
        return 0;
    return 1;
}

uint64_t Size(void *matrix)
{
    return CsxSize<index_t, value_t>(matrix);
}

uint64_t SizeSym(void *matrix)
{
    return CsxSymSize<index_t, value_t>(matrix);
}

void DestroyCsx(void *matrix)
{
    spm_mt_t *spm_mt = (spm_mt_t *)matrix;
    PutSpmMt<index_t, value_t>(spm_mt);
}

}
