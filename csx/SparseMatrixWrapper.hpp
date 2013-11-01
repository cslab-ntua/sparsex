/*
 * SparseMatrixWrapper.hpp -- Wrappers of the SparseMatrix routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSE_MATRIX_WRAPPER_HPP
#define SPARSE_MATRIX_WRAPPER_HPP

#ifdef __cplusplus

#include "../C-API/mattype.h"

extern "C" {
#endif

/**
 *  SparseMatrix wrapper functions.
 */
void *CreateCSR(index_t *rowptr, index_t *colind, value_t *values,
                index_t nr_rows, index_t nr_cols, int zero_based);
void *CreateMMF(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz);
void *ReorderCSR(void *matrix, index_t **permutation);
void *ReorderMMF(void *matrix, index_t **permutation);
void *TuneCSR(void *matrix, int *symmetric);
void *TuneMMF(void *matrix, int *symmetric);
void DestroyCSR(void *matrix);
void DestroyMMF(void *matrix);
void SaveTuned(void *matrix, const char *filename, index_t *permutation);
void *LoadTuned(const char *filename, index_t *nr_rows, index_t *nr_cols,
                index_t *nnz, index_t **permutation);
int GetValue(void *matrix, index_t row, index_t col, value_t *value);
int SetValue(void *matrix, index_t row, index_t col, value_t value);
/* uint64_t Size(void *matrix); */
/* uint64_t SizeSym(void *matrix); */
void DestroyCsx(void *matrix);
void SetPropertyByMnemonic(const char *key, const char *value);
void SetPropertiesFromEnv();

#ifdef __cplusplus
}
#endif

#endif // SPARSE_MATRIX_WRAPPER_HPP
