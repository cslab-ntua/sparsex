/**
 * \file mat_vec.h -- \brief Sparse matrix routines.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef SPARSEX_MAT_VEC_H
#define SPARSEX_MAT_VEC_H

#include "common.h"

/**
 *  Creates and returns a valid tunable matrix object from a Compressed 
 *  Sparse Row (CSR) representation.
 *
 *  @param[in] rowptr       array \a rowptr of the CSR format.
 *  @param[in] colind       array \a colind of the CSR format.
 *  @param[in] values       array \a values of the CSR format.
 *  @param[in] nr_rows      number of rows of the matrix.
 *  @param[in] nr_cols      number of columns of the matrix.
 *  @param[in]              optional argument that specifies the indexing
 *                          (either \c INDEXING_ZERO_BASED or 
 *                          \c INDEXING_ONE_BASED).
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
 *  This routine releases any memory internally used by the sparse input
 *  handle \a input.
 *
 *  @param[in] input        the input matrix handle.
 *  @return                 an error code.
 */
spx_error_t spx_input_destroy(spx_input_t *input);

/**
 *  \brief Converts the input matrix into the CSX format.
 *
 *  @param[in] input        the input matrix handle.
 *  @param[in] ...          optional flag that indicates whether the matrix 
 *                          should be reordered with use of the Reverse Cuthill
 *                          McKee algorithm (\c OP_REORDER).
 *  @return                 a handle to the tuned matrix.
 */
spx_matrix_t *spx_mat_tune(spx_input_t *input, ...);

/**
 *  Returns the value of the corresponding element in (\a row, \a column),
 *  where \a row and \a column can be either zero- or one-based indexes. If the
 *  element doesn't exist an error is returned.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] row          the a row of the element to be retrieved.
 *  @param[in] column       the column of the element to be retrieved.
 *  @param[out] value       the value of the element in (\a row, \a column).
 *  @param[in]              optional argument that specifies the indexing
 *                          (either \c INDEXING_ZERO_BASED or 
 *                          \c INDEXING_ONE_BASED).
 *  @return                 an error code.
 */
spx_error_t spx_mat_get_entry(const spx_matrix_t *A, spx_index_t row,
                              spx_index_t column, spx_value_t *value, ...);

/**
 *  Sets the value of the corresponding element in (\a row, \a column), where
 *  \a row and \a column can be either zero- or one-based indexes. If the element
 *  doesn't exist an error is returned.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @param[in] row          the row of the element to be set.
 *  @param[in] column       the column of the element to be set.
 *  @param[in] value        the new value of the element in (\a row, \a column).
 *  @param[in]              optional argument that specifies the indexing
 *                          (either \c INDEXING_ZERO_BASED or 
 *                          \c INDEXING_ONE_BASED).
 *  @return                 an error code.
 */
spx_error_t spx_mat_set_entry(spx_matrix_t *A, spx_index_t row,
                              spx_index_t column, spx_value_t value, ...);

/**
 *  \brief Stores the matrix in the CSX format into a binary file.
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
 *  \brief Returns the number of rows of the matrix.
 */
spx_index_t spx_mat_get_nrows(const spx_matrix_t *A);

/**
 *  \brief Returns the number of columns of the matrix.
 */
spx_index_t spx_mat_get_ncols(const spx_matrix_t *A);

/**
 *  \brief Returns the number of non-zeros of the matrix.
 */
spx_index_t spx_mat_get_nnz(const spx_matrix_t *A);

/**
 *  \brief Returns information on the partitioning of the matrix.
 */
spx_partition_t *spx_mat_get_parts(spx_matrix_t *A);

/**
 *  \brief Returns the permutation applied to the matrix.
 */
spx_perm_t *spx_mat_get_perm(const spx_matrix_t *A);

/**
 *  This routines performs a matrix-vector multiplication of the following form:
 *                      <em> y <-- alpha*A*x + beta*y </em>
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
spx_error_t spx_matvec_mult(spx_scalar_t alpha, const spx_matrix_t *A,
                            spx_vector_t *x, spx_scalar_t beta, spx_vector_t *y);

/**
 *  This routines performs a matrix-vector multiplication of the following form:
 *                      <em> y <-- alpha*A*x + beta*y </em>
 *  where \a alpha and \a beta are scalars, \a x and \a y are vectors
 *  and \a A is a sparse matrix. The matrix is originally given in the CSR
 *  format and converted internally into the CSX format. This higher-level routine
 *  hides the preprocessing phase of CSX.
 *
 *  This routine can be efficiently used in a loop, since only the first call
 *  will convert the matrix into the CSX format and every subsequent call will use
 *  the previously tuned matrix handle.
 *
 *  @param[in] A            either an invalid matrix handle or a tuned
 *                          matrix handle. If A is equal to an \c INVALID_MAT then
 *                          the matrix in the CSR format is first converted to CSX. 
 *                          Otherwise, the valid (previously) tuned matrix handle
 *                          is used to perform the multiplication. 
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
spx_error_t spx_csr_matvec_mult(spx_matrix_t *A, 
                                spx_index_t nr_rows, spx_index_t nr_cols,
                                spx_index_t *rowptr, spx_index_t *colind, 
                                spx_value_t *values,
                                spx_scalar_t alpha, spx_vector_t *x, 
                                spx_scalar_t beta, spx_vector_t *y);

/**
 *  This routine releases any memory internally used by the tuned matrix
 *  handle \a A.
 *
 *  @param[in] A            the tuned matrix handle.
 *  @return                 an error code.
 */
spx_error_t spx_mat_destroy(spx_matrix_t *A);

/**
 *  This routine creates a partitioning handle of the matrix in the Compressed
 *  Sparse Row (CSR) format.
 *
 *  @param[in] rowptr       array \a rowptr of the CSR format.
 *  @param[in] nr_rows      number of rows of the matrix.
 *  @param[in] nr_threads   number of partitions of the matrix.
 *  @return                 an error code.
 */
spx_partition_t *spx_partition_csr(spx_index_t *rowptr, spx_index_t nr_rows, 
                                   unsigned int nr_threads);

/**
 *  This routine releases any memory internally used by the partition
 *  handle \a p.
 *
 *  @param[in] p            the partitioning handle.
 *  @return                 an error code.
 */
spx_error_t spx_partition_destroy(spx_partition_t *p);

/**
 *  Sets the \a option according to the description in \a string for the tuning
 *  process to follow. For available tuning options \see common.h
 *
 *  @param[in] option       the option to be set.
 *  @param[in] string       a description of how to set the option.
 */
void spx_option_set(const char *option, const char *string);

/**
 *  Sets the tuning options according to the environmental variables 
 *  set on the command line. For available tuning options \see common.h
 */
void spx_options_set_from_env();

/**
 *  Creates and returns a valid vector object, whose values must be explicitly
 *  initialized.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @param[in] p            a partitioning handle.
 *  @return                 a valid vector object.
 */
spx_vector_t *vec_create(unsigned long size, void *p);
spx_vector_t *vec_create_numa(unsigned long size, spx_partition_t *p);
#ifdef SPM_NUMA
#   define spx_vec_create vec_create_numa
#else
#   define spx_vec_create vec_create
#endif

/**
 *  Creates and returns a valid vector object, whose values are mapped to a
 *  user-defined array.
 *
 *  @param[in] buff         the user-supplied buffer.
 *  @param[in] size         the size of the buffer.
 *  @param[in] p            a partitioning handle.
 *  @return                 a valid vector object.
 */
spx_vector_t *vec_create_from_buff(spx_value_t *buff, unsigned long size,
                                   void *p);
spx_vector_t *vec_create_from_buff_numa(spx_value_t *buff, unsigned long size,
                                        spx_partition_t *p);
#ifdef SPM_NUMA
#   define spx_vec_create_from_buff vec_create_from_buff_numa
#else
#   define spx_vec_create_from_buff vec_create_from_buff
#endif

/**
 *  Creates and returns a valid vector object, whose values are randomly filled.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @param[in] p            a partitioning handle.
 *  @return                 a valid vector object.
 */
spx_vector_t *vec_create_random(unsigned long size, void *p);
spx_vector_t *vec_create_random_numa(unsigned long size, spx_partition_t *p);
#ifdef SPM_NUMA
#   define spx_vec_create_random vec_create_random_numa
#else
#   define spx_vec_create_random vec_create_random
#endif

/**
 *  Initializes the valid vector object \a v with \a val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 */
void vec_init(spx_vector_t *v, spx_value_t val);
#define spx_vec_init vec_init

/**
 *  Initializes the [\a start, \a end) part of the valid vector object \a v
 *  with \a val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_init_part(spx_vector_t *v, spx_value_t val, spx_index_t start,
                   spx_index_t end);
#define spx_vec_init_part vec_init_part

/**
 *  Initializes the valid vector object \a v with random values in the
 *  range [\a min, \a max].
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] max          maximum value of initializing range.
 *  @param[in] min          minimum value of initializing range.
 */
void vec_init_rand_range(spx_vector_t *v, spx_value_t max, spx_value_t min);
#define spx_vec_init_rand_range vec_init_rand_range

/**
 *  Sets the element at index \a idx of vector \a v to be equal to \a val.
 *  \a idx is assumed to be one-based.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] idx          an index inside the vector.
 *  @param[in] val          the value to be set.
 */
spx_error_t spx_vec_set_entry(spx_vector_t *v, spx_index_t idx, spx_value_t val);

/**
 *  \brief v2 -> num * v1
 *
 *  Scales the input vector \a v1 by a constant value \a num and places the 
 *  result in vector \a v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] num          the const by which to scale \a v1.
 */
void vec_scale(spx_vector_t *v1, spx_vector_t *v2, spx_scalar_t num);
#define spx_vec_scale vec_scale

/**
 *  \brief v3 -> v1 + num * v2
 *
 *  Scales the input vector \a v2 by a constant value \a num, adds the
 *  result to vector \a v1 and places the result in vector \a v3.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] num          the scalar by which to scale \a v1.
 */
void vec_scale_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3, 
                   spx_scalar_t num);
#define spx_vec_scale_add vec_scale_add

/**
 *  \brief \a v3[\a start...\a end-1] -> \a v1[\a start...\a end-1] + \a num * \a v2[\a start...\a end-1]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] num          the scalar by which to scale \a v1.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_scale_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                        spx_scalar_t num, spx_index_t start, spx_index_t end);
#define spx_vec_scale_add_part vec_scale_add_part

/**
 *  \brief v3 -> v1 + v2
 *
 *  Adds the input vectors \a v1 and \a v2 and places the result in \a v3.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 */
void vec_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3);
#define spx_vec_add vec_add

/**
 *  \brief v3[start...end-1] -> v1[start...end-1] + v2[start...end-1]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                  spx_index_t start, spx_index_t end);
#define spx_vec_add_part vec_add_part

/**
 *  \brief v3 -> v1 - v2
 *
 *  Subtracts the input vector \a v2 from \a v1 and places the result in \a v3.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 */
void vec_sub(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3);
#define spx_vec_sub vec_sub

/**
 *  \brief v3[start...end-1] -> v1[start...end-1] - v2[start...end-1]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_sub_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                  spx_index_t start, spx_index_t end);
#define spx_vec_sub_part vec_sub_part

/**
 *  Returns the product of the input vectors \a v1 and \a v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @return                 the product of the input vectors.
 */
spx_value_t vec_mul(const spx_vector_t *v1, const spx_vector_t *v2);
#define spx_vec_mul vec_mul

/**
 *  \brief Returns v1[start...end-1] * v2[start...end-1]
 *
 *  Returns the product of the part [\a start, \a end) of the input vectors
 *  \a v1 and \a v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 *  @return                 the product of the input vectors.
 */
spx_value_t vec_mul_part(const spx_vector_t *v1, const spx_vector_t *v2,
                         spx_index_t start, spx_index_t end);
#define spx_vec_mul_part vec_mul_part

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
 *  Copies the elements of v1 to v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
void vec_copy(const spx_vector_t *v1, spx_vector_t *v2);
#define spx_vec_copy vec_copy

/**
 *  Compares the elements of v1 and v2. If they are equal it returns 0,
 *  else -1.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
int vec_compare(const spx_vector_t *v1, const spx_vector_t *v2);
#define spx_vec_compare vec_compare

/**
 *  \brief Prints the input vector v.
 *
 *  @param[in] v            a valid vector object.
 *  @return                 error code.
 */
void vec_print(const spx_vector_t *v);
#define spx_vec_print vec_print

/**
 *  \brief Destroys the input vector v.
 *
 *  @param[in] v            a valid vector object.
 */
void vec_destroy(spx_vector_t *v);
#define spx_vec_destroy vec_destroy


spx_error_t spx_matvec_mult_pool(spx_scalar_t alpha, const spx_matrix_t *A, 
                                 spx_vector_t *x, spx_scalar_t beta, spx_vector_t *y);
spx_thread_pool_t *spx_thread_pool_create(size_t nthreads);

#endif // SPARSEX_MAT_VEC_H

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
