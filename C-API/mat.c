/**
 * libcsx/mat.c -- Sparse matrix routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "mat.h"
#include "CsxMatvec.hpp"
#include "SparseMatrixWrapper.hpp"

#include <unistd.h>
#include <float.h>

#define SPMV_FN matvec_mt
#define SPMV_SYM_FN matvec_sym_mt

/**
 *  \brief The sparse matrix internal representation.
 */
struct matrix {
    index_t nrows, ncols, nnz;
    perm_t *permutation;    /**< The permutation, in case the matrix has been
                               reordered */
    void *csx;              /**< The tuned matrix representation, i.e. the input
                               matrix transformed to the CSX format */
    int symmetric;
};

/**
 *  \brief The input matrix internal representation.
 */
struct input {
    index_t nrows, ncols, nnz;
    char type;   /**< A character that stores the format of the input matrix
                    ('C' for CSR or 'M' for MMF) */
    void *mat;   /**< The input matrix representation */
};

/**
 *  \brief Allocates a matrix structure.
 */
static void mat_alloc_struct(matrix_t **A)
{
    *A = libcsx_malloc(matrix_t, sizeof(matrix_t));
    (*A)->nrows = 0;
    (*A)->ncols = 0;
    (*A)->nnz = 0;
    (*A)->permutation = INVALID_PERM;
    (*A)->csx = INVALID_MAT;
    (*A)->symmetric = 0;
}

/**
 *  \brief Frees a matrix structure.
 */
static void mat_free_struct(matrix_t *A)
{
    if (A->csx) {
        DestroyCsx(A->csx);
    }
    if (A->permutation) {
        libcsx_free(A->permutation);
    }
    libcsx_free(A);
}

/**
 *  \brief Allocates an input structure.
 */
static void input_alloc_struct(input_t **A)
{
    *A = libcsx_malloc(input_t, sizeof(input_t));
    (*A)->nrows = 0;
    (*A)->ncols = 0;
    (*A)->nnz = 0;
    (*A)->type = 'U';
    (*A)->mat = INVALID_INPUT;
}

/**
 *  \brief Frees an input structure.
 */
static void input_free_struct(input_t *A)
{
    if (A->mat) {
        if (A->type == 'C') {
            DestroyCSR(A->mat);
        } else {
            DestroyMMF(A->mat);
        }
    }
    libcsx_free(A);
}

input_t *libcsx_mat_create_csr(index_t *rowptr, index_t *colind, value_t *values,
                               index_t nr_rows, index_t nr_cols, int zero_based)
{
    /* Check validity of input arguments */
    if (!check_dim(nr_rows) || !check_dim(nr_cols)) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix dimensions");
        return INVALID_INPUT;
    }
    if (!check_indexing(zero_based)) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid indexing");
        return INVALID_INPUT;
    }
    if (!rowptr) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid rowptr argument");
        return INVALID_INPUT;
    }
    if (!colind) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid colind argument");
        return INVALID_INPUT;
    }
    if (!values) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid values argument");
        return INVALID_INPUT;
    }

    /* Create CSR wrapper */
    input_t *A = INVALID_INPUT;
    input_alloc_struct(&A);
    A->type = 'C';
    A->nrows = nr_rows;
    A->ncols = nr_cols;
    A->nnz = rowptr[nr_rows] - !zero_based;
    A->mat = CreateCSR(rowptr, colind, values, nr_rows, nr_cols, zero_based);
    if (!A->mat) {
        SETERROR_1(LIBCSX_ERR_INPUT_MAT, "creating CSR wrapper failed");
        return INVALID_INPUT;
    }

    return A;
}

input_t *libcsx_mat_create_mmf(const char *filename)
{
    /* Check validity of input argument */
    if (!filename) {
        SETERROR_0(LIBCSX_ERR_FILE);
        return INVALID_INPUT;
    }
    if (access(filename, F_OK | R_OK) == -1) {  //not 100% safe check
        SETERROR_0(LIBCSX_ERR_FILE);
        return INVALID_INPUT;
    }

    /* Load matrix */
    input_t *A = INVALID_INPUT;
    input_alloc_struct(&A);
    A->type = 'M';
    A->mat = CreateMMF(filename, &(A->nrows), &(A->ncols), &(A->nnz));
    if (!A->mat) {
        SETERROR_1(LIBCSX_ERR_INPUT_MAT, "loading matrix from MMF file failed");
        return INVALID_INPUT;
    }

    return A;
}

matrix_t *libcsx_mat_tune(input_t *input, int reorder)
{
    /* Check validity of input argument */
    if (!input) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid input matrix");
        return INVALID_MAT;
    }

    /* Convert to CSX */
    matrix_t *tuned = INVALID_MAT;
    perm_t *permutation = INVALID_PERM;
    mat_alloc_struct(&tuned);
    tuned->nrows = input->nrows;
    tuned->ncols = input->ncols;
    tuned->nnz = input->nnz;
    if (input->type == 'C') {
        if (reorder) {
            input->mat = ReorderCSR(input->mat, &permutation);
        }
        tuned->csx = TuneCSR(input->mat, &(tuned->symmetric));
    } else if (input->type == 'M') {
        if (reorder) {
            input->mat = ReorderMMF(input->mat, &permutation);
        }
        tuned->csx = TuneMMF(input->mat, &(tuned->symmetric));
    }
    tuned->permutation = permutation;
    permutation = INVALID_PERM;

    return tuned;
}

libcsx_error_t libcsx_matvec_mult(const matrix_t *A, value_t alpha, vector_t *x,
                                  value_t beta, vector_t *y)
{
    /* Check validity of input arguments */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_ARG_INVALID;
    }
    if (!x) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid vector x");
        return LIBCSX_ERR_ARG_INVALID;
    }
    if (!y) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid vector y");
        return LIBCSX_ERR_ARG_INVALID;
    }

    /* Check compatibility between matrix and vector dimensions */
    if (!check_vec_dim(x, A->ncols) && !check_vec_dim(y, A->nrows)) {
        SETERROR_0(LIBCSX_ERR_VEC_DIM);
        return LIBCSX_ERR_VEC_DIM;
    }

    /* Reorder vectors */
    perm_t *p = A->permutation;
    vector_t *permuted_x = x;
    vector_t *permuted_y = y;
    if (p != INVALID_PERM) {
        permuted_x = vec_reorder(x, p);
        permuted_y = vec_reorder(y, p);
    }
    /* Compute kernel */
    if (!A->symmetric) {
        SPMV_FN(A->csx, permuted_x, alpha, permuted_y, beta);
    } else {
        SPMV_SYM_FN(A->csx, permuted_x, alpha, permuted_y, beta);
    }
    /* Restore resulting vector to original order */
    if (p != INVALID_PERM) {
        vector_t *tmp_y = vec_inv_reorder(permuted_y, p);
        vec_copy(tmp_y, y);
        vec_destroy(tmp_y);
    }

    /* Cleanup */
    if (p != INVALID_PERM) {
        vec_destroy(permuted_x);
        permuted_x = INVALID_VEC;
        vec_destroy(permuted_y);
        permuted_y = INVALID_VEC;
    }
    p = INVALID_PERM;

    return LIBCSX_SUCCESS;
}

libcsx_error_t libcsx_mat_get_entry(const matrix_t *A, index_t row,
                                    index_t column, value_t *value)
{
    /* Check validity of input argument */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_TUNED_MAT;
    }

    /* Search for element */
    if (A->permutation != INVALID_PERM) {
        row = A->permutation[row - 1] + 1;
        column = A->permutation[column - 1] + 1;
    }
    int err = GetValue(A->csx, row, column, value);
    if (err != 0) {
        SETERROR_0(LIBCSX_ERR_ENTRY_NOT_FOUND);
        return LIBCSX_ERR_ENTRY_NOT_FOUND;
    }

    return LIBCSX_SUCCESS;
}

libcsx_error_t libcsx_mat_set_entry(matrix_t *A, index_t row, index_t column,
                                    value_t value)
{
    /* Check validity of input argument */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_TUNED_MAT;
    }

    /* Set new value */
    if (A->permutation != INVALID_PERM) {
        row = A->permutation[row - 1] + 1;
        column = A->permutation[column - 1] + 1;
    }
    int err = SetValue(A->csx, row, column, value);
    if (err != 0) {
        SETERROR_0(LIBCSX_ERR_ENTRY_NOT_FOUND);
        return LIBCSX_ERR_ENTRY_NOT_FOUND;
    }

    return LIBCSX_SUCCESS;
}

libcsx_error_t libcsx_mat_save(const matrix_t *A, const char *filename)
{
    /* Check validity of input arguments */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_TUNED_MAT;
    }
    if (!filename) {
        SETWARNING(LIBCSX_WARN_CSXFILE);
    }

    /* Save the tuned matrix */
    SaveTuned(A->csx, filename, A->permutation);

    return LIBCSX_SUCCESS;
}

matrix_t *libcsx_mat_restore(const char *filename)
{
    /* Check validity of input argument */
    if (!filename) {
        SETERROR_0(LIBCSX_ERR_FILE);
        return INVALID_MAT;
    }
    if (access(filename, F_OK | R_OK) == -1) {  //not 100% safe check
        SETERROR_0(LIBCSX_ERR_FILE);
        return INVALID_MAT;
    }

    /* Load tuned matrix from file */
    matrix_t *A = INVALID_MAT;
    mat_alloc_struct(&A);
    A->csx = LoadTuned(filename, &(A->nrows), &(A->ncols), &(A->nnz),
                       &(A->permutation));

    return A;
}

index_t libcsx_mat_get_nrows(const matrix_t *A)
{
    return  A->nrows;
}

index_t libcsx_mat_get_ncols(const matrix_t *A)
{
    return A->ncols;
}

libcsx_error_t libcsx_mat_destroy_tuned(matrix_t *A)
{
    /* Check validity of input argument */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_TUNED_MAT;
    }

    /* Free allocated memory of matrix handle */
    mat_free_struct(A);

    return LIBCSX_SUCCESS;
}

libcsx_error_t libcsx_mat_destroy_input(input_t *A)
{
    /* Check validity of input argument */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid input handle");
        return LIBCSX_ERR_INPUT_MAT;
    }

    /* Free allocated memory of matrix handle */
    input_free_struct(A);

    return LIBCSX_SUCCESS;
}

void libcsx_set_option(const char *option, const char *value)
{
    SetPropertyByMnemonic(option, value);
}
