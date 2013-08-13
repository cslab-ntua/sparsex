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
#include "SparseMatrixWrapper.h"
#include <unistd.h>
#include <float.h>

#ifdef SPM_NUMA
    #include "spmv_loops_mt_numa.h"
    #include "spmv_loops_sym_mt_numa.h"
    #define SPMV_FN spmv_double_bench_mt_loop_numa
    #define SPMV_SYM_FN spmv_double_bench_sym_mt_loop_numa
#else
    #include "spmv_loops_mt.h"
// #include "../lib/prfcnt/tsc.h"
//#include "spmv_loops_sym_mt.h"
    #define SPMV_FN spmv_double_mt
    #define SPMV_SYM_FN spmv_double_sym_mt
#endif // SPM_NUMA

static double calc_imbalance(void *m);
static matrix_t *tune(input_t *input);
static vector_t *mult(const matrix_t *A, const vector_t *x);

/**
 *  \brief The sparse matrix internal representation.
 */
struct matrix {
    index_t nrows, ncols, nnz;
    int symmetric;
    perm_t *permutation;   /**< The permutation, in case the matrix has been reordered */
    void *tuned;    /**< The tuned matrix representation, i.e. the input matrix tuned in the CSX format */
    double pre_time;
};

/**
 *  \brief The input matrix internal representation.
 */
struct input {
    index_t nrows, ncols, nnz;
    unsigned int tuning_options;
    void *mat;   /**< The input matrix representation */
    char type;   /**< A character that stores the format of the input matrix ('C' for CSR or 'M' for MMF) */
};

/**
 *  \brief Allocates a matrix structure.
 */
static libcsx_error_t mat_alloc_struct(matrix_t **A)
{
    *A = libcsx_malloc(matrix_t, sizeof(matrix_t));
    (*A)->nrows = 0;
    (*A)->ncols = 0;
    (*A)->nnz = 0;
    (*A)->tuned = NULL;
    (*A)->permutation = NULL;
    (*A)->symmetric = 0;
    (*A)->pre_time = 0;

    return LIBCSX_SUCCESS;
}

/**
 *  \brief Frees a matrix structure.
 */
static libcsx_error_t mat_free_struct(matrix_t *A)
{
    if (A->tuned) {
        DestroyCsx(A->tuned);
    }
    if (A->permutation) {
        libcsx_free(A->permutation);
    }
    libcsx_free(A);

    return LIBCSX_SUCCESS;
}

/**
 *  \brief Allocates an input structure.
 */
static libcsx_error_t input_alloc_struct(input_t **A)
{
    *A = libcsx_malloc(input_t, sizeof(input_t));
    (*A)->type = 'U';
    (*A)->nrows = 0;
    (*A)->ncols = 0;
    (*A)->nnz = 0;
    (*A)->mat = NULL;
    (*A)->tuning_options = 0;

    return LIBCSX_SUCCESS;
}

/**
 *  \brief Frees an input structure.
 */
libcsx_error_t input_free_struct(input_t *A)
{
    if (A->mat) {
        if (A->type == 'C') {
            DestroyCSR(A->mat);
        } else {
            DestroyMMF(A->mat);
        }
    }
    libcsx_free(A);

    return LIBCSX_SUCCESS;
}

input_t *mat_create_csr(index_t *rowptr, index_t *colind, value_t *values,
                        index_t nr_rows, index_t nr_cols, int zero_based)
{
    /* Check validity of input arguments */
    if (!check_dim(nr_rows) || !check_dim(nr_cols)) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix dimensions");
        return INVALID_MAT;
    }
    if (!check_indexing(zero_based)) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid indexing");
        return INVALID_MAT;
    }
    if (!rowptr) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid rowptr argument");
        return INVALID_MAT;
    }
    if (!colind) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "Invalid colind argument");
        return INVALID_MAT;
    }
    if (!values) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid values argument");
        return INVALID_MAT;
    }

    /* Create CSR wrapper */
    input_t *A = INVALID_MAT;
    input_alloc_struct(&A);
    A->type = 'C';
    A->nrows = nr_rows;
    A->ncols = nr_cols;
    A->nnz = rowptr[nr_rows] - !zero_based;
    A->mat = CreateCSR(rowptr, colind, values, nr_rows, nr_cols, zero_based);
    if (!A->mat) {
        SETERROR_1(LIBCSX_ERR_INPUT_MAT, "creating CSR wrapper failed");
        return INVALID_MAT;
    }

    return A;
}

input_t *mat_create_mmf(const char *filename)
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

    /* Load matrix */
    input_t *A = INVALID_MAT;
    input_alloc_struct(&A);
    A->type = 'M';
    A->mat = CreateMMF(filename, &(A->nrows), &(A->ncols), &(A->nnz));
    if (!A->mat) {
        SETERROR_1(LIBCSX_ERR_INPUT_MAT, "loading matrix from MMF file failed");
        return INVALID_MAT;
    }

    return A;
}

matrix_t *mat_tune(input_t *input)
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

vector_t *matvec_mult(const matrix_t *A, value_t alpha, const vector_t *x,
                      value_t beta, vector_t *y)
{
    /* Check validity of input arguments */
    if (!x) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid vector x");
        return INVALID_VEC;
    }

    /* Check compatibility between matrix and vector dimensions */
    if (!check_vec_dim(x, A->ncols)) {
        SETERROR_0(LIBCSX_ERR_VEC_DIM);
        return INVALID_VEC;
    }
    if (y && !check_vec_dim(y, A->nrows)) {
        SETERROR_0(LIBCSX_ERR_VEC_DIM);
        return INVALID_VEC;
    }

    /* x <-- permuted(x) */
    perm_t *p = A->permutation;
    vector_t *permuted_x = INVALID_VEC;
    vector_t *permuted_y = y;
    if (p != NULL) {
        permuted_x = vec_reorder(x, p);
        if (y != NULL) {
            permuted_y = vec_reorder(y, p);
            vec_destroy(y);
            y = INVALID_VEC;
        }
    }
    /* x <-- alpha*x */
    vector_t *scaled_x = (vector_t *) x;
    if (permuted_x != NULL) {
        scaled_x = permuted_x;
    }
    if (alpha != 1) {
        if (permuted_x != NULL) {
            scaled_x = vec_scale_original(permuted_x, alpha);
        } else {
            scaled_x = vec_scale(x, alpha);
        }
    }
    /* y <-- beta*y */
    vector_t *scaled_y = permuted_y;
    if (beta != 1 && permuted_y != NULL) {
        scaled_y = vec_scale_original(permuted_y, beta);
    }
    /* mv_result <-- alpha*A*x */
    vector_t *mv_result = INVALID_VEC;
    mv_result = mult(A, scaled_x);
    /* y <-- alpha*A*x + beta*y*/
    if (scaled_y != NULL) {
        mv_result = vec_add(mv_result, scaled_y);
    }
    /* Restore to original order */
    if (p != NULL) {
        y = vec_inv_reorder(mv_result, p);
        vec_destroy(mv_result);
    } else {
        y = mv_result;
    }

    /* Cleanup */
    if (permuted_x != NULL) {
        vec_destroy(permuted_x);
        permuted_x = INVALID_VEC;
        vec_destroy(permuted_y);
        permuted_y = INVALID_VEC;
    } else if (alpha != 1) {
        vec_destroy(scaled_x);
    }
    scaled_x = INVALID_VEC;
    scaled_y = INVALID_VEC;
    mv_result = INVALID_VEC;
    p = INVALID_PERM;

    return y;
}

value_t mat_get_entry(const matrix_t *A, index_t row, index_t column)
{
    /* Check validity of input argument */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_TUNED_MAT;
    }

    /* Search for element */
    if (A->permutation != NULL) {
        row = A->permutation[row - 1] + 1;
        column = A->permutation[column - 1] + 1;
    }
    value_t value = GetValue(A->tuned, row, column);
    if (value == 0) {                                   //this isn't correct
        SETERROR_0(LIBCSX_ERR_ENTRY_NOT_FOUND);
        return LIBCSX_ERR_ENTRY_NOT_FOUND;
    }

    return value;
}

libcsx_error_t mat_set_entry(matrix_t *A, index_t row, index_t column,
                             value_t value)
{
    /* Check validity of input argument */
    if (!A) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid matrix handle");
        return LIBCSX_ERR_TUNED_MAT;
    }

    /* Set new value */
    if (A->permutation != NULL) {
        row = A->permutation[row - 1] + 1;
        column = A->permutation[column - 1] + 1;
    }
    if (!SetValue(A->tuned, row, column, value)) {              //check this
        SETERROR_0(LIBCSX_ERR_ENTRY_NOT_FOUND);
        return LIBCSX_ERR_ENTRY_NOT_FOUND;
    }

    return LIBCSX_SUCCESS;
}

libcsx_error_t mat_save(const matrix_t *A, const char *filename)
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
    if (A->permutation != NULL) {
        SaveTuned(A->tuned, filename, A->permutation);
    } else {
        SaveTuned(A->tuned, filename, NULL);
    }

    return LIBCSX_SUCCESS;
}

matrix_t *mat_restore(const char *filename)
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

index_t mat_get_nrows(const matrix_t *A)
{
    return  A->nrows;
}

index_t mat_get_ncols(const matrix_t *A)
{
    return A->ncols;
}

libcsx_error_t mat_destroy_tuned(matrix_t *A)
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

libcsx_error_t mat_destroy_input(input_t *A)
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

// OPT_SYMMETRIC and OPT_REORDER can only be set (not cleared) for the time being
void set_tuning_option(input_t *input, option_t option, const char *string)
{
    switch(option) {
    case OPT_SYMMETRIC:
        input->tuning_options |= OPT_SYMMETRIC;        
        break;
    case OPT_REORDER:
        input->tuning_options |= OPT_REORDER;        
        break;
    case OPT_SPLIT_BLOCKS:
//    set_property(?, string);
        break;
    case OPT_ONE_DIM_BLOCKS:
//    set_property(?, string);
        break;
    case OPT_XFORMS:
//    set_property(?, string);
        break;
    case OPT_WINDOW_SIZE:
//    set_property(?, string);
        break;
    case OPT_NR_SAMPLES:
//    set_property(?, string);
        break;
    case OPT_SAMPLING_PORTION:
//    set_property(?, string);
        break;
    default:
        SETERROR_0(LIBCSX_WARN_TUNING_OPT);
    }
}

void set_runtime_option(option_t option, const char *string)
{
    switch(option) {
    case OPT_NR_THREADS:
//    set_property(?, string);
    case OPT_AFFINITY:
//    set_property(?, string);
    default:
        SETERROR_0(LIBCSX_WARN_RUNTIME_OPT);
    }
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
//    else error

    return A;
}

static vector_t *mult(const matrix_t *A, const vector_t *x)
{
    VECTOR_TYPE *y = NULL;

    y = VECTOR_NAME(_create)(A->nrows);
    if (!A->symmetric) {
        double secs = SPMV_FN(A->tuned, (VECTOR_TYPE *) x, y, 1);
        double flops = (double) (A->nnz*2) / ((double) 1000*1000*secs);
        printf("m:%s s:%lu pt:%lf t:%lf r:%lf\n", "csx", Size(A->tuned),
               A->pre_time, secs, flops);
        double imbalance = calc_imbalance(A->tuned);
        printf("Load imbalance: %lf%s\n", 100*imbalance, "%");
    } else {
        /* double secs = SPMV_SYM_FN(A->tuned, (VECTOR_TYPE *) x, y); */
        /* double flops = (double) (A->nnz*2) / ((double) 1000*1000*secs); */
        /* printf("m:%s t:%lf r:%lf\n", "csx-sym", secs, flops); */
        /* printf("m:%s ms:%lu s:%lu pt:%lf t:%lf r:%lf\n", "csx-sym",  */
        /*        map_size(A->tuned), SizeSym(A->tuned), pre_time, secs, */
        /*        flops); */
    }

    return (vector_t *) y;
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
