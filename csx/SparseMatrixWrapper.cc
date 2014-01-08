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

#include "CsxGetSet.hpp"
#include "CsxSaveRestore.hpp"
#include "CsxUtil.hpp"
#include "Runtime.hpp"
#include "SparseMatrix.hpp"
#include "SparseMatrixWrapper.hpp"

/**
 *  Explicit instantiation definitions.
 */
template class SparseMatrix<CSR<index_t, value_t> >;
template class SparseMatrix<MMF<index_t, value_t> >;

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
    return (void *) matrix;
}

void *CreateMMF(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz)
{
    SparseMatrix<MMF<index_t, value_t> > *matrix =
        new SparseMatrix<MMF<index_t, value_t> >(filename);
    *nr_rows = static_cast<index_t>(matrix->GetNrRows());
    *nr_cols = static_cast<index_t>(matrix->GetNrCols());
    *nnz = static_cast<index_t>(matrix->GetNrNonzeros());
    
    return (void *) matrix;
}

void *ReorderCSR(void *matrix, index_t **permutation)
{
    SparseMatrix<CSR<index_t, value_t> > *mat = 
        (SparseMatrix<CSR<index_t, value_t> > *) matrix;
    std::vector<size_t> perm;

    mat->Reorder(perm);
    if (!perm.empty()) {
        *permutation = new index_t[mat->GetNrRows()];
        std::copy(perm.begin(), perm.end(), *permutation);
    }

    return (void *) mat;
}

void *ReorderMMF(void *matrix, index_t **permutation) 
{
    SparseMatrix<MMF<index_t, value_t> > *mat = 
        (SparseMatrix<MMF<index_t, value_t> > *) matrix;
    std::vector<size_t> perm;

    mat->Reorder(perm);
    if (!perm.empty()) {
        *permutation = new index_t[mat->GetNrRows()];
        std::copy(perm.begin(), perm.end(), *permutation);
    }

    return (void *) mat;
}

void *TuneCSR(void *matrix, int *symmetric)
{
    SparseMatrix<CSR<index_t, value_t> > *mat = 
        (SparseMatrix<CSR<index_t, value_t> > *) matrix;
    RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    rt_context.SetRuntimeContext(config);

    spm_mt_t *spm_mt = mat->CreateCsx();
    if (config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric)) {
        *symmetric = 1;
    }

    return (void *) spm_mt;
}

void *TuneMMF(void *matrix, int *symmetric)
{
    SparseMatrix<MMF<index_t, value_t> > *mat = 
        (SparseMatrix<MMF<index_t, value_t> > *) matrix;
    RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();
    RuntimeContext &rt_context = RuntimeContext::GetInstance();
    rt_context.SetRuntimeContext(config);

    spm_mt_t *spm_mt = mat->CreateCsx();
    if (config.GetProperty<bool>(RuntimeConfiguration::MatrixSymmetric)) {
        *symmetric = 1;
    }

    return (void *) spm_mt;
}

void DestroyCSR(void *matrix)
{
    SparseMatrix<CSR<index_t, value_t> > *mat = 
        (SparseMatrix<CSR<index_t, value_t> > *) matrix;
    delete mat;
}

void DestroyMMF(void *matrix)
{
    SparseMatrix<MMF<index_t, value_t> > *mat = 
        (SparseMatrix<MMF<index_t, value_t> > *) matrix;
    delete mat;
}

void SaveTuned(void *matrix, const char *filename, index_t *permutation)
{
    spm_mt_t *spm_mt = (spm_mt_t *) matrix;
    SaveCsx<index_t, value_t>(spm_mt, filename, permutation);
}

void *LoadTuned(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz, index_t **permutation)
{
    spm_mt_t *spm_mt = 0;
    csx_t<value_t> *csx = 0;

    spm_mt = RestoreCsx<index_t, value_t>(filename, permutation);
	for (unsigned int i = 0; i < spm_mt->nr_threads; i++) {
        csx = (csx_t<value_t> *) spm_mt->spm_threads[i].spm;
        *nr_rows += csx->nrows;
        *nnz += csx->nnz;
    }
    *nr_cols = csx->ncols;

    return (void *) spm_mt;
}

int GetValue(void *matrix, index_t row, index_t col, value_t *value)
{
    spm_mt_t *spm_mt = (spm_mt_t *) matrix;
    if (GetValueCsx<index_t, value_t>(spm_mt, row, col, value))
        return 0;
    return -1;
}

int SetValue(void *matrix, index_t row, index_t col, value_t value)
{
    spm_mt_t *spm_mt = (spm_mt_t *) matrix;
    if (SetValueCsx<index_t, value_t>(spm_mt, row, col, value))
        return 0;
    return -1;
}

void DestroyCsx(void *matrix)
{
    spm_mt_t *spm_mt = (spm_mt_t *) matrix;
    PutSpmMt<value_t>(spm_mt);
}

void SetPropertyByMnemonic(const char *key, const char *value)
{
    RuntimeConfiguration &rt_config = RuntimeConfiguration::GetInstance();
    rt_config.SetProperty(rt_config.GetMnemonic(string(key)), string(value));
}

void SetPropertiesFromEnv()
{
    RuntimeConfiguration &config = RuntimeConfiguration::GetInstance();
    config.LoadFromEnv();
}

// uint64_t Size(void *matrix)
// {
//     return CsxSize<value_t>(matrix);
// }

// uint64_t SizeSym(void *matrix)
// {
//     return CsxSymSize<value_t>(matrix);
// }

}
