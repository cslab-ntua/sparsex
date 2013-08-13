/*
 * SparseMatrixWrapper.h -- Wrappers of the SparseMatrix routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSEMATRIXWRAPPER_H__
#define SPARSEMATRIXWRAPPER_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "../C-API/mattype.h"

/**
 *  SparseMatrix wrapper functions.
 */
void *CreateCSR(index_t *rowptr, index_t *colind, value_t *values,
                index_t nr_rows, index_t nr_cols, int zero_based);
void *CreateMMF(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz);
void *ReorderCSR(void *matrix, index_t **permutation);
void *ReorderMMF(void *matrix, index_t **permutation);
void *TuneCSR(void *matrix, double *time);
void *TuneMMF(void *matrix, double *time);
void DestroyCSR(void *matrix);
void DestroyMMF(void *matrix);
void SaveTuned(void *matrix, const char *filename, index_t *permutation);
void *LoadTuned(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz, int *symmetric, index_t **permutation,
                double *pre_time);
value_t GetValue(void *matrix, index_t row, index_t col);
int SetValue(void *matrix, index_t row, index_t col, value_t value);
uint64_t Size(void *matrix);
uint64_t SizeSym(void *matrix);
void DestroyCsx(void *matrix);

#ifdef __cplusplus
}
#endif

#endif // SPARSEMATRIXWRAPPER_H__
