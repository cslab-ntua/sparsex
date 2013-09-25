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
#include "csx_matvec.h"
#include "SparseMatrixWrapper.h"
#include <unistd.h>
#include <float.h>

#define SPMV_FN matvec_mt
#define SPMV_SYM_FN matvec_sym_mt

static double calc_imbalance(void *m);
static matrix_t *tune(input_t *input);

/**
 *  \brief The sparse matrix internal representation.
 */
struct matrix {
    index_t nrows, ncols, nnz;
    int symmetric;
    perm_t *permutation;    /**< The permutation, in case the matrix has been
                               reordered */
    void *tuned;            /**< The tuned matrix representation, i.e. the input
                               matrix tuned in the CSX format */
    double pre_time;
};

/**
 *  \brief The input matrix internal representation.
 */
struct input {
    index_t nrows, ncols, nnz;
    unsigned int tuning_options;
    void *mat;   /**< The input matrix representation */
    char type;   /**< A character that stores the format of the input matrix
                    ('C' for CSR or 'M' for MMF) */
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
    (*A)->tuned = INVALID_MAT;
    (*A)->permutation = INVALID_PERM;
    (*A)->symmetric = 0;
    (*A)->pre_time = 0;
}

/**
 *  \brief Frees a matrix structure.
 */
static void mat_free_struct(matrix_t *A)
{
    if (A->tuned) {
        DestroyCsx(A->tuned);
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
    (*A)->type = 'U';
    (*A)->nrows = 0;
    (*A)->ncols = 0;
    (*A)->nnz = 0;
    (*A)->mat = INVALID_INPUT;
    (*A)->tuning_options = 0;
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

matrix_t *libcsx_mat_tune(input_t *input)
{
    /* Check validity of input argument */
    if (!input) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid input matrix");
        return INVALID_MAT;
    }

    /* Convert to CSX */
    matrix_t *tuned = INVALID_MAT;
    perm_t *permutation = INVALID_PERM;
    if (input->tuning_options & OPT_REORDER) {
        if (input->type == 'C') {
            input->mat = ReorderCSR(input->mat, &permutation);
        } else if (input->type == 'M') {
            input->mat = ReorderMMF(input->mat, &permutation);
        }
    }
    tuned = tune(input);
    tuned->permutation = permutation;
    permutation = INVALID_PERM;
    if (!tuned) {
        SETERROR_0(LIBCSX_ERR_TUNED_MAT);
        return INVALID_MAT;
    }

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
        SPMV_FN(A->tuned, permuted_x, alpha, permuted_y, beta);
    } else {
        /* SPMV_SYM_FN(A->tuned, permuted_x, alpha, permuted_y, beta); */
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
    int err = GetValue(A->tuned, row, column, value);
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
    int err = SetValue(A->tuned, row, column, value);
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
    SaveTuned(A->tuned, filename, A->permutation);

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
    A->tuned = LoadTuned(filename, &(A->nrows), &(A->ncols), &(A->nnz),
                         &(A->symmetric), &(A->permutation), &(A->pre_time));

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

libcsx_error_t libcsx_set_tuning_option(input_t *input, option_t option,
                                        const char *string)
{
    switch(option) {
    case OPT_SYMMETRIC:
        input->tuning_options |= OPT_SYMMETRIC;        
        break;
    case OPT_REORDER:
        input->tuning_options |= OPT_REORDER;        
        break;
    case OPT_SPLIT_BLOCKS:
        break;
    case OPT_ONE_DIM_BLOCKS:
        break;
    case OPT_XFORMS:
        break;
    case OPT_WINDOW_SIZE:
        break;
    case OPT_NR_SAMPLES:
        break;
    case OPT_SAMPLING_PORTION:
        break;
    default:
        SETERROR_0(LIBCSX_WARN_TUNING_OPT);
    }

    return LIBCSX_SUCCESS;
}

libcsx_error_t libcsx_set_runtime_option(option_t option, const char *string)
{
    switch(option) {
    case OPT_NR_THREADS:
    case OPT_AFFINITY:
    default:
        SETERROR_0(LIBCSX_WARN_RUNTIME_OPT);
    }

    return LIBCSX_SUCCESS;
}

static matrix_t *tune(input_t *input)
{
    matrix_t *A = NULL;

    mat_alloc_struct(&A);
    A->nrows = input->nrows;
    A->ncols = input->ncols;
    A->nnz = input->nnz;
    if (input->tuning_options & OPT_SYMMETRIC) {
        A->symmetric = 1;
    }
    if (input->type == 'C') {
        A->tuned = TuneCSR(input->mat, &(A->pre_time));
    } else if (input->type == 'M') {
        A->tuned = TuneMMF(input->mat, &(A->pre_time));
    }

    return A;
}

static double calc_imbalance(void *A)
{
	spm_mt_t *spm_mt = (spm_mt_t *) A;
	size_t i;

	double min_time = DBL_MAX;
	double max_time = 0.0;
	double total_time = 0.0;
	size_t worst = -1;
	for (i = 0; i < spm_mt->nr_threads; ++i) {
		spm_mt_thread_t *spm = &(spm_mt->spm_threads[i]);
		double thread_time = spm->secs;
		printf("Thread %zd time: %lf\n", i, thread_time);
		total_time += thread_time;
		if (thread_time > max_time) {
			max_time = thread_time;
			worst = i;
		}

		if (thread_time < min_time)
			min_time = thread_time;
	}

	double ideal_time = total_time / spm_mt->nr_threads;
	printf("Worst thread: %zd\n", worst);
	printf("Expected perf. improvement: %.2f %%\n",
	       100*(max_time / ideal_time - 1));
	return (max_time - min_time) / min_time;
}
