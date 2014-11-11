/*
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013-2014, Athena Elafrou
 * Copyright (C) 2014,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file matvec.h
 * \brief Sparse matrix routines
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#ifndef SPARSEX_MATVEC_H
#define SPARSEX_MATVEC_H

#include <sparsex/common.h>

/**
 *  Creates and returns a valid tunable matrix object from a Compressed 
 *  Sparse Row (CSR) representation.
 *
 *  @param[in] rowptr       array \a rowptr of the CSR format.
 *  @param[in] colind       array \a colind of the CSR format.
 *  @param[in] values       array \a values of the CSR format.
 *  @param[in] nr_rows      number of rows of the matrix.
 *  @param[in] nr_cols      number of columns of the matrix.
 *  @param[in] ...          optional argument that specifies the indexing
 *                          (either \c SPX_INDEX_ZERO_BASED or 
 *                          \c SPX_INDEX_ONE_BASED).
 *  @return                 a handle to the input matrix.
 */
spx_input_t *spx_input_load_csr(spx_index_t *rowptr, spx_index_t *colind,
                                spx_value_t *values, spx_index_t nr_rows, 
                                spx_index_t nr_cols, ...);

/**
 *  Creates and returns a valid tunable matrix object from a file in the
 *  Matrix Market File Format (MMF).
 *
 *  @param[in] filename     name of the file where the matrix is stored.
 *  @return                 a handle to the input matrix.
 */
spx_input_t *spx_input_load_mmf(const char *filename);

/**
 *  Releases any memory internally used by the sparse input
 *  handle \a input.
 *
 *  @param[in] input        the input matrix handle.
 *  @return                 an error code.
 */
spx_error_t spx_input_destroy(spx_input_t *input);

/**
 *  Converts the input matrix into the CSX format by applying all the options
 *  previously set with the spx_option_set() routine. In case no options
 *  have been explicitly set, the default values are used (see Table 3.2 of the
 *  User's Guide).
 *
 *  @param[in] input        the input matrix handle.
 *  @param[in] ...          optional flag that indicates whether the matrix 
 *                          should be reordered with use of the Reverse Cuthill
 *                          McKee algorithm (\c SPX_MAT_REORDER).
 *  @return                 a handle to the tuned matrix.
 */
spx_matrix_t *spx_mat_tune(spx_input_t *input, ...);

/**
 *  Returns the value of the corresponding nonzero element in
 *  (\a row, \a column), where \a row and \a column can be either zero-
 *  or one-based indices. Default indexing is zero-based, but it can be
 *  overridden through the optional flag. If the element exists, its value
 *  is returned in \a value.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] row          the a row of the element to be retrieved.
 *  @param[in] column       the column of the element to be retrieved.
 *  @param[out] value       the value of the element in (\a row, \a column).
 *  @param[in] ...          optional argument that specifies the indexing
 *                          (either \c SPX_INDEX_ZERO_BASED or 
 *                          \c SPX_INDEX_ONE_BASED).
 *  @return                 an error code.
 */
spx_error_t spx_mat_get_entry(const spx_matrix_t *A, spx_index_t row,
                              spx_index_t column, spx_value_t *value, ...);

/**
 *  Sets the value of the corresponding element in (\a row, \a column), where \a
 *  row and \a column can be either zero- or one-based indices. If the element
 *  doesn't exist an error is returned.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] row          the row of the element to be set.
 *  @param[in] column       the column of the element to be set.
 *  @param[in] value        the new value of the element in (\a row, \a column).
 *  @param[in] ...          optional argument that specifies the indexing
 *                          (either \c SPX_INDEX_ZERO_BASED or 
 *                          \c SPX_INDEX_ONE_BASED).
 *  @return                 an error code.
 */
spx_error_t spx_mat_set_entry(spx_matrix_t *A, spx_index_t row,
                              spx_index_t column, spx_value_t value, ...);

/**
 *  Stores the matrix in the CSX format into a binary file.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] filename     the name of the binary file where the matrix
 *                          will be dumped.
 *  @return                 an error code.
 */
spx_error_t spx_mat_save(const spx_matrix_t *A, const char *filename);

/**
 *  Reconstructs the matrix in the CSX format from a binary file.
 *
 *  @param[in] filename     the name of the file where the matrix is stored.
 *  @return                 a handle to the tuned matrix.
 */
spx_matrix_t *spx_mat_restore(const char *filename);

/**
 *  Returns the number of rows of the matrix.
 */
spx_index_t spx_mat_get_nrows(const spx_matrix_t *A);

/**
 *  Returns the number of columns of the matrix.
 */
spx_index_t spx_mat_get_ncols(const spx_matrix_t *A);

/**
 *  Returns the number of non-zeros of the matrix.
 */
spx_index_t spx_mat_get_nnz(const spx_matrix_t *A);

/**
 *  Returns a partitioning object for matrix \a A.
 */
spx_partition_t *spx_mat_get_partition(spx_matrix_t *A);

/**
 *  Returns the permutation computed for matrix \a A by applying the Reverse
 *  Cuthill-McKee algorithm.
 */
spx_perm_t *spx_mat_get_perm(const spx_matrix_t *A);

/**
 *  Performs a matrix-vector multiplication of the following form:
 *                      <em> y = alpha*A*x </em>
 *  where \a alpha is a scalar, \a x and \a y are vectors
 *  and \a A is a sparse matrix in the CSX format.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] alpha        a scalar.
 *  @param[in] x            the input vector.
 *  @param[in,out] y        the output vector.
 *  @return                 an error code.
 */
spx_error_t spx_matvec_mult(spx_value_t alpha, const spx_matrix_t *A,
                            spx_vector_t *x, spx_vector_t *y);

/**
 *  Performs a matrix-vector multiplication of the following form:
 *                      <em> y = alpha*A*x + beta*y </em>
 *  where \a alpha and \a beta are scalars, \a x and \a y are vectors
 *  and \a A is a sparse matrix in the CSX format.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] alpha        a scalar.
 *  @param[in] x            the input vector.
 *  @param[in] beta         a scalar.
 *  @param[in,out] y        the output vector.
 *  @return                 an error code.
 */
spx_error_t spx_matvec_kernel(spx_value_t alpha, const spx_matrix_t *A,
                              spx_vector_t *x, spx_value_t beta,
                              spx_vector_t *y);

/**
 *  Performs a matrix-vector multiplication of the following form:
 *                      <em> y = alpha*A*x + beta*y </em>
 *  where \a alpha and \a beta are scalars, \a x and \a y are vectors
 *  and \a A is a sparse matrix. The matrix is originally given in the CSR
 *  format and converted internally into the CSX format. This higher-level
 *  routine hides the preprocessing phase of CSX.
 *
 *  This routine can be efficiently used in a loop, since only the first call
 *  will convert the matrix into the CSX format and every subsequent call will
 *  use the previously tuned matrix handle.
 *
 *  @param[in] A            either a pointer to an invalid matrix handle or a
 *                          tuned matrix handle. If (*A) is equal to
 *                          \c SPX_INVALID_MAT then the matrix in the CSR
 *                          format is first converted to CSX. Otherwise,
 *                          the valid (previously) tuned matrix handle is
 *                          used to perform the multiplication.
 *  @param[in] nr_rows      number of rows of the matrix \a A.
 *  @param[in] nr_cols      number of columns of the matrix \a A.
 *  @param[in] rowptr       array \a rowptr of the CSR format.
 *  @param[in] colind       array \a colind of the CSR format.
 *  @param[in] values       array \a values of the CSR format.
 *  @param[in] alpha        a scalar.
 *  @param[in] x            the input vector.
 *  @param[in] beta         a scalar.
 *  @param[in,out] y        the output vector.
 *  @return                 an error code.
 */
spx_error_t spx_matvec_kernel_csr(spx_matrix_t **A, 
                                  spx_index_t nr_rows, spx_index_t nr_cols,
                                  spx_index_t *rowptr, spx_index_t *colind, 
                                  spx_value_t *values,
                                  spx_value_t alpha, spx_vector_t *x, 
                                  spx_value_t beta, spx_vector_t *y);

/**
 *  Releases any memory internally used by the tuned matrix
 *  handle \a A.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @return                 an error code.
 */
spx_error_t spx_mat_destroy(spx_matrix_t *A);

/**
 *  Creates a partitioning object for the matrix in the Compressed
 *  Sparse Row (CSR) format. This routine should be used in conjunction
 *  with the spx_matvec_kernel_csr() multiplication routine.
 *
 *  @param[in] rowptr       array \a rowptr of the CSR format.
 *  @param[in] nr_rows      number of rows of the matrix.
 *  @param[in] nr_threads   number of partitions of the matrix.
 *  @return                 an error code.
 */
spx_partition_t *spx_partition_csr(spx_index_t *rowptr, spx_index_t nr_rows, 
                                   size_t nr_threads);

/**
 *  Releases any memory internally used by the partition
 *  handle \a p.
 *
 *  @param[in] p            the partitioning handle.
 *  @return                 an error code.
 */
spx_error_t spx_partition_destroy(spx_partition_t *p);

/**
 *  Sets the \a option according to the description in \a string for the tuning
 *  process to follow. For available tuning options and how to set them refer
 *  to Table 3.2 of the User's Guide.
 *
 *  @param[in] option       the option to be set.
 *  @param[in] string       a description of how to set the option.
 */
void spx_option_set(const char *option, const char *string);

/**
 *  Sets the tuning options according to the environmental variables 
 *  set on the command line (to be deprecated).
 */
void spx_options_set_from_env();

/**
 *  Creates and returns a valid vector object, whose values are initialized to
 *  zero.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @param[in] p            a partitioning handle.
 *  @return                 a valid vector object.
 */
spx_vector_t *spx_vec_create(size_t size, spx_partition_t *p);

/**
 *  Creates and returns a valid vector object, whose values are mapped to a
 *  user-defined array. If SPX_VEC_AS_IS is set, then the input buffer will
 *  be shared with the library and modifications will directly apply to it.
 *  In this case pointers \a buff and \a tuned point to the same memory
 *  location. If SPX_VEC_TUNE is selected, the buffer provided by the user
 *  might be copied into an optimally allocated buffer (depending on the
 *  platform) and \a tuned might point to this buffer. Thus, the user must
 *  always check whether \a buff equals \a tuned. If the buffer is actually
 *  tuned, then it should be used instead of the original. The common free()
 *  function applies to this buffer and will have to be explicitly called by
 *  the user.
 *
 *  @param[in] buff         the user-supplied buffer.
 *  @param[in] tuned        the tuned buffer.
 *  @param[in] size         the size of the buffer.
 *  @param[in] p            a partitioning handle.
 *  @param[in] mode         the vector mode (either \c SPX_VEC_AS_IS or
 *                          \c SPX_VEC_TUNE).
 *  @return                 a valid vector object.
 */
spx_vector_t *spx_vec_create_from_buff(spx_value_t *buff,
                                       spx_value_t **tuned,
                                       size_t size, spx_partition_t *p,
                                       spx_vecmode_t mode);

/**
 *  Creates and returns a valid vector object, whose values are randomly filled.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @param[in] p            a partitioning handle.
 *  @return                 a valid vector object.
 */
spx_vector_t *spx_vec_create_random(size_t size, spx_partition_t *p);

/**
 *  Initializes the vector object \a v with \a val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 */
void spx_vec_init(spx_vector_t *v, spx_value_t val);

/**
 *  Initializes the [\a start : \a end] range of the vector object \a v
 *  with \a val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void spx_vec_init_part(spx_vector_t *v, spx_value_t val, spx_index_t start,
                       spx_index_t end);

/**
 *  Initializes the vector object \a v with random values in the
 *  range [\a min, \a max].
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] max          maximum value of initializing range.
 *  @param[in] min          minimum value of initializing range.
 */
void spx_vec_init_rand_range(spx_vector_t *v, spx_value_t max, spx_value_t min);

/**
 *  Sets the element at index \a idx of vector \a v to be equal to \a val.
 *  \a idx is assumed to be one-based.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] idx          an index inside the vector.
 *  @param[in] val          the value to be set.
 *  @param[in] ...          optional argument that specifies the indexing
 *                          (either \c SPX_INDEX_ZERO_BASED or 
 *                          \c SPX_INDEX_ONE_BASED).

 */
spx_error_t spx_vec_set_entry(spx_vector_t *v, spx_index_t idx,
                              spx_value_t val, ...);

/**
 *
 *  Scales the input vector \a v1 by a constant value \a num and places the 
 *  result in vector \a v2.
 *  \a v2 -> \a num * \a v1
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] num          the constant by which to scale \a v1.
 */
void spx_vec_scale(spx_vector_t *v1, spx_vector_t *v2, spx_value_t num);

/**
 *  Scales the input vector \a v2 by a constant value \a num, adds the
 *  result to vector \a v1 and places the result in vector \a v3.
 *  \a v3 -> \a v1 + \a num * \a v2
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] num          the scalar by which to scale \a v1.
 */
void spx_vec_scale_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3, 
                       spx_value_t num);

/**
 *  \a v3[\a start : \a end] = \a v1[\a start : \a end] +
 *  \a num * \a v2[\a start : \a end]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] num          the scalar by which to scale \a v1.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void spx_vec_scale_add_part(spx_vector_t *v1, spx_vector_t *v2,
                            spx_vector_t *v3, spx_value_t num,
                            spx_index_t start, spx_index_t end);

/**
 *  Adds the input vectors \a v1 and \a v2 and places the result in \a v3.
 *  \a v3 = \a v1 + \a v2
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 */
void spx_vec_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3);

/**
 *  Adds the range [\a start : \a end] of the input vectors \a v1 and \a v2 and
 *  places the result in \a v3.
 *  \a v3[\a start : \a end] = \a v1[\a start : \a end] +
 *  \a v2[\a start : \a end]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void spx_vec_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                      spx_index_t start, spx_index_t end);

/**
 *  Subtracts the input vector \a v2 from \a v1 and places the result in \a v3.
 *  \a v3 -> \a v1 - \a v2
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 */
void spx_vec_sub(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3);

/**
 *  Subtracts the input vector \a v2 from \a v1 in the range [\a start :
 *  \a end] and places the result in \a v3.
 *  \a v3[\a start : \a end] = \a v1[\a start : \a end] -
 *  \a v2[\a start : \a end]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void spx_vec_sub_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                      spx_index_t start, spx_index_t end);

/**
 *  Returns the product of the input vectors \a v1 and \a v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @return                 the product of the input vectors.
 */
spx_value_t spx_vec_mul(const spx_vector_t *v1, const spx_vector_t *v2);

/**
 *  Returns the product of the range [\a start : \a end] of the input vectors
 *  \a v1 and \a v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 *  @return                 the product of the input vectors.
 */
spx_value_t spx_vec_mul_part(const spx_vector_t *v1, const spx_vector_t *v2,
                             spx_index_t start, spx_index_t end);

/**
 *  Reorders the input vector v according to the permutation p, leaving
 *  the original vector intact.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] p            a permutation.
 *  @return                 the permuted input vector.
 */
spx_error_t spx_vec_reorder(spx_vector_t *v, spx_perm_t *p);

/**
 *  Inverse-reorders the permuted input vector \a v1, according
 *  to the permutation \a p, storing the result in \a v2.
 *
 *  @param[in,out] v        the vector object to be inverse-reordered.
 *  @param[in] p            a permutation.
 *  @return                 an error code.
 */
spx_error_t spx_vec_inv_reorder(spx_vector_t *v, spx_perm_t *p);

/**
 *  Copies the values of v1 to v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
void spx_vec_copy(const spx_vector_t *v1, spx_vector_t *v2);

/**
 *  Compares the values of v1 and v2. If they are equal it returns 0,
 *  else -1.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @return                 0 on success, -1 on failure.
 */
int spx_vec_compare(const spx_vector_t *v1, const spx_vector_t *v2);

/**
 *  Prints the input vector v.
 *
 *  @param[in] v            a valid vector object.
 */
void spx_vec_print(const spx_vector_t *v);

/**
 *  Destroys the input vector v, after optionally moving data to a user-
 *  defined buffer.
 *
 *  @param[in] v            a valid vector object.
 */
void spx_vec_destroy(spx_vector_t *v);

#endif  /* SPARSEX_MATVEC_H */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
