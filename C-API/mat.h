/**
 * libcsx/mat.h -- Sparse matrix routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_MAT_H__
#define LIBCSX_MAT_H__

#include "common.h"
#include "vec.h"

/**
 *  \brief The sparse matrix handle type.
 */
typedef struct matrix matrix_t;

/**
 *  \brief The input matrix handle type.
 */
typedef struct input input_t;

/**
 *  Creates and returns a valid tunable matrix object from a Compressed 
 *  Sparse Row (CSR) representation.
 *
 *  @param[in] rowptr       array "rowptr" of CSR format.
 *  @param[in] colind       array "colind" of CSR format.
 *  @param[in] values       array "values" of CSR format.
 *  @param[in] nr_rosw      number of rows.
 *  @param[in] nr_cols      number of columns.
 *  @param[in] zero_based   indexing.
 *  @return                 a handle to the input matrix.
 */
input_t *mat_create_csr(index_t *rowptr, index_t *colind, value_t *values,
                        index_t nr_rows, index_t nr_cols, int zero_based);

/**
 *  Creates and returns a valid tunable matrix object from a file in the
 *  Matrix Market Format (MMF).
 *
 *  @param[in] filename     name of the file where the matrix is kept.
 *  @return                 a handle to the input matrix.
 */
input_t *mat_create_mmf(const char *filename);

/**
 *  Converts the input matrix in the CSX format according to the options set.
 *
 *  @param[in] input        the input matrix.
 *  @return                 a handle to the tuned matrix.
 */
matrix_t *mat_tune(input_t *input);

/* Implements y <-- alpha*A*x + beta*y */
vector_t *matvec_mult(const matrix_t *A, value_t alpha, const vector_t *x,
                      value_t beta, vector_t *y);

/**
 *  Returns the value of the corresponding element in (row, column), where
 *  row and column are assumed to be one-based indexes. If the element doesn't
 *  exist an error is returned.
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] row          the row of the element to be retrieved.
 *  @param[in] column       the column of the element to be retrieved.
 *  @return                 the value of the element (row, column).
 */
value_t mat_get_entry(const matrix_t *A, index_t row, index_t column);

/**
 *  Sets the value of the corresponding element in (row, column), where
 *  row and column are assumed to be one-based indexes. If the element doesn't
 *  exist an error is returned.
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] row          the row of the element to be set.
 *  @param[in] column       the column of the element to be set.
 *  @param[out] value       the new value of the element in (row, column).
 *  @return                 an error code.
 */
libcsx_error_t mat_set_entry(matrix_t *A, index_t row, index_t column,
                             value_t value);

/**
 *  Dumps the matrix in the CSX format in a binary file.
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] filename     name of the file where the matrix will be dumped.
 *  @return                 an error code.
 */
libcsx_error_t mat_save(const matrix_t *A, const char *filename);

/**
 *  Reconstructs the matrix in the CSX format from a binary file.
 *
 *  @param[in] filename     name of the file where the matrix is dumped.
 *  @return                 a handle to the tuned matrix.
 */
matrix_t *mat_restore(const char *filename);

/**
 *  This routine #FIXME Utility for drawing sparse matrices.
 *  
 *
 *  @param[in] A            the tuned matrix.
 *  @return                 an error code.
 */
libcsx_error_t mat_draw(const matrix_t *A);

/**
 *  \brief Returns the number of rows of the matrix.
 */
index_t mat_get_nrows(const matrix_t *A);

/**
 *  \brief Returns the number of columns of the matrix.
 */
index_t mat_get_ncols(const matrix_t *A);

/**
 *  This routine releases any memory internally used by the sparse matrix
 *  handle A.
 *
 *  @param[in] A            the tuned matrix.
 *  @return                 an error code.
 */
libcsx_error_t mat_destroy_tuned(matrix_t *A);

/**
 *  This routine releases any memory internally used by the sparse input
 *  handle A.
 *
 *  @param[in] A            the input matrix.
 *  @return                 an error code.
 */
libcsx_error_t mat_destroy_input(input_t *A);

/**
 *  Sets the option #option according to the string #string for the tuning
 *  process to follow. For available tuning options \see #FIXME libcsx/common.h.
 *
 *  @param[in] A            the input matrix to be tuned.
 *  @param[in] option       the option to be set.
 *  @param[in] string       a description of how to set the option.
 */
void set_tuning_option(input_t *A, option_t option , const char *string);

/**
 *  Sets the option #option according to the string #string.
 *  For available runtime options \see #FIXME libcsx/common.h
 *
 *  @param[in] option       the option to be set.
 *  @param[in] string       a description of how to set the option.
 */
void set_runtime_option(option_t option , const char *string);

#endif // LIBCSX_MAT_H__

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
