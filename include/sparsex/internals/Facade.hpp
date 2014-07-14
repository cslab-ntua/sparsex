/*
 * \file Facade.hpp
 *
 * \brief Wrappers of the SparseMatrix routines
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#ifndef SPARSEX_INTERNALS_FACADE_HPP
#define SPARSEX_INTERNALS_FACADE_HPP

#include <sparsex/internals/cdecl.h>
#include <sparsex/types.h>

SPX_BEGIN_C_DECLS__

/**
 *  SparseMatrix wrapper functions.
 */
void *CreateCSR(spx_index_t *rowptr, spx_index_t *colind, spx_value_t *values,
                spx_index_t nr_rows, spx_index_t nr_cols, int zero_based);
void *CreateMMF(const char *filename,
                spx_index_t *nr_rows, spx_index_t *nr_cols, spx_index_t *nnz);
void *ReorderCSR(void *matrix, spx_index_t **permutation);
void *ReorderMMF(void *matrix, spx_index_t **permutation);
void *TuneCSR(void *matrix, int *symmetric);
void *TuneMMF(void *matrix, int *symmetric);
void DestroyCSR(void *matrix);
void DestroyMMF(void *matrix);
void SaveTuned(void *matrix, const char *filename, spx_index_t *permutation);
void *LoadTuned(const char *filename,
                spx_index_t *nr_rows, spx_index_t *nr_cols,
                spx_index_t *nnz, spx_index_t **permutation);
int GetValue(void *matrix,
             spx_index_t row, spx_index_t col, spx_value_t *value);
int SetValue(void *matrix, spx_index_t row, spx_index_t col, spx_value_t value);
void DestroyCsx(void *matrix);
void SetPropertyByMnemonic(const char *key, const char *value);
void SetPropertiesFromEnv();
void GetNodes(int *nodes);
void CreatePool();

SPX_END_C_DECLS__

#endif // SPARSEX_INTERNALS_FACADE_HPP
