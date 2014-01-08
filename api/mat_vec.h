/**
 * SparseX/mat_vec.h -- Sparse matrix routines.
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
 *  @param[in] rowptr       array "rowptr" of CSR format.
 *  @param[in] colind       array "colind" of CSR format.
 *  @param[in] values       array "values" of CSR format.
 *  @param[in] nr_rows      number of rows.
 *  @param[in] nr_cols      number of columns.
 *  @param[in] indexing     indexing.
 *  @return                 a handle to the input matrix.
 */
input_t *spx_input_load_csr(index_t *rowptr, index_t *colind, value_t *values,
                            index_t nr_rows, index_t nr_cols, 
                            property_t indexing);

/**
 *  Creates and returns a valid tunable matrix object from a file in the
 *  Matrix Market File Format (MMF).
 *
 *  @param[in] filename     name of the file where the matrix is kept.
 *  @return                 a handle to the input matrix.
 */
input_t *spx_input_load_mmf(const char *filename);

/**
 *  This routine releases any memory internally used by the sparse input
 *  handle A.
 *
 *  @param[in] A            the input matrix.
 *  @return                 an error code.
 */
spx_error_t spx_input_destroy(input_t *A);

/**
 *  Converts the input matrix into the CSX format.
 *
 *  @param[in] input        the input matrix.
 *  @param[in] ...          optional flag that indicates whether the matrix 
 *                          should be reordered, <OP_REORDER>.
 *  @return                 a handle to the tuned matrix.
 */
matrix_t *spx_mat_tune(input_t *input, ...);

/**
 *  Returns the value of the corresponding element in (row, column), where
 *  row and column are assumed to be one-based indexes. If the element doesn't
 *  exist an error is returned.
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] row          the row of the element to be retrieved.
 *  @param[in] column       the column of the element to be retrieved.
 *  @param[out] value       the value of the element in (row, column).
 *  @return                 an error code.
 */
spx_error_t spx_mat_get_entry(const matrix_t *A, index_t row,
                              index_t column, value_t *value);

/**
 *  Sets the value of the corresponding element in (row, column), where
 *  row and column are assumed to be one-based indexes. If the element doesn't
 *  exist an error is returned.
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] row          the row of the element to be set.
 *  @param[in] column       the column of the element to be set.
 *  @param[in] value        the new value of the element in (row, column).
 *  @return                 an error code.
 */
spx_error_t spx_mat_set_entry(matrix_t *A, index_t row, index_t column,
                              value_t value);

/**
 *  Stores the matrix in the CSX format into a binary file.
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] filename     the name of the file where the matrix will be dumped.
 *  @return                 an error code.
 */
spx_error_t spx_mat_save(const matrix_t *A, const char *filename);

/**
 *  Reconstructs the matrix in the CSX format from a binary file.
 *
 *  @param[in] filename     the name of the file where the matrix is dumped.
 *  @return                 a handle to the tuned matrix.
 */
matrix_t *spx_mat_restore(const char *filename);

/**
 *  \brief Returns the number of rows of the matrix.
 */
index_t spx_mat_get_nrows(const matrix_t *A);

/**
 *  \brief Returns the number of columns of the matrix.
 */
index_t spx_mat_get_ncols(const matrix_t *A);

/**
 *  \brief Returns information on the partitioning of the matrix.
 */
partition_t *spx_mat_get_parts(matrix_t *A);

/**
 *  \brief Returns the permutation applied to the matrix.
 */
perm_t *spx_mat_get_perm(const matrix_t *A);

/**
 *  This routines performs a matrix-vector multiplication of the following form:
 *                      y <-- alpha*A*x + beta*y
 *
 *  @param[in] A            the tuned matrix.
 *  @param[in] alpha        a scalar.
 *  @param[in] x            the input vector.
 *  @param[in] beta         a scalar.
 *  @param[in/out] y        the output vector.
 *  @return                 an error code.
 */
spx_error_t spx_matvec_mult(scalar_t alpha, const matrix_t *A, vector_t *x,
                            scalar_t beta, vector_t *y);

/**
 *  This routine releases any memory internally used by the sparse matrix
 *  handle A.
 *
 *  @param[in] A            the tuned matrix.
 *  @return                 an error code.
 */
spx_error_t spx_mat_destroy(matrix_t *A);

/**
 *  This routine releases any memory internally used by the partition
 *  handle p.
 *
 *  @param[in] p            the partitioning handle.
 *  @return                 an error code.
 */
spx_error_t spx_part_destroy(partition_t *p);

/**
 *  Sets the option #option according to the string #string for the tuning
 *  process to follow. For available tuning options \see #FIXME
 *
 *  @param[in] A            the input matrix to be tuned.
 *  @param[in] option       the option to be set.
 *  @param[in] string       a description of how to set the option.
 */
void spx_option_set(const char *option, const char *string);

/**
 *  Sets the tuning options according to the environmental variables 
 *  set on the command line. For available tuning options \see #FIXME
 */
void spx_options_set_from_env();

/**
 *  Creates and returns a valid vector object, whose values must be explicitly
 *  initialized.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @return                 a valid vector object.
 */
vector_t *vec_create(unsigned long size, void *p);
vector_t *vec_create_numa(unsigned long size, partition_t *p);
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
 *  @return                 a valid vector object.
 */
vector_t *vec_create_from_buff(value_t *buff, unsigned long size, void *p);
vector_t *vec_create_from_buff_numa(double *buff, unsigned long size, 
                                    partition_t *p);
#ifdef SPM_NUMA
#   define spx_vec_create_from_buff vec_create_from_buff_numa
#else
#   define spx_vec_create_from_buff vec_create_from_buff
#endif

/**
 *  Creates and returns a valid vector object, whose values are randomly filled.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @return                 a valid vector object.
 */
vector_t *vec_create_random(unsigned long size, void *p);
vector_t *vec_create_random_numa(unsigned long size, partition_t *p);
#ifdef SPM_NUMA
#   define spx_vec_create_random vec_create_random_numa
#else
#   define spx_vec_create_random vec_create_random
#endif

/**
 *  Initializes the valid vector object @v with @val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 */
void vec_init(vector_t *v, double val);
#define spx_vec_init vec_init

/**
 *  Initializes the [@start, @end) part of the valid vector object @v
 *  with @val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void spxvec_init_part(vector_t *v, double val, unsigned long start,
                         unsigned long end);
#define spx_vec_init_part vec_init_part

/**
 *  Initializes the valid vector object @v with random values in the
 *  range [@min, @max].
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] max          maximum value of initializing range.
 *  @param[in] min          minimum value of initializing range.
 */
void vec_init_rand_range(vector_t *v, double max, double min);
#define spx_vec_init_rand_range vec_init_rand_range

/**
 *  Sets the element at index @idx of vector @v to be equal to @val.
 *  @idx is assumed to be one-based.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] idx          an index inside the vector.
 *  @param[in] val          the value to be set.
 */
spx_error_t spx_vec_set_entry(vector_t *v, int idx, double val);

/**
 *  \brief v2 -> num * v1
 *
 *  Scales the input vector @v1 by a constant value @num and places the 
 *  result in vector @v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] num          the const by which to scale @v1.
 */
void vec_scale(vector_t *v1, vector_t *v2, double num);
#define spx_vec_scale vec_scale

/**
 *  \brief v3 -> v1 + num * v2
 *
 *  Scales the input vector @v2 by a constant value @num, adds the
 *  result to vector @v1 and places the result in vector @v3.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] num          the const by which to scale @v1..
 */
void vec_scale_add(vector_t *v1, vector_t *v2, vector_t *v3, double num);
#define spx_vec_scale_add vec_scale_add

/**
 *  \brief v3[start...end-1] -> v1[start...end-1] + num * v2[start...end-1]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_scale_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                        double num, unsigned long start,
                        unsigned long end);
#define spx_vec_scale_add_part vec_scale_add_part

/**
 *  \brief v3 -> v1 + v2
 *
 *  Adds the input vectors @v1 and @v2 and places the result in @v3.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 */
void vec_add(vector_t *v1, vector_t *v2, vector_t *v3);
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
void vec_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                  unsigned long start, unsigned long end);
#define spx_vec_add_part vec_add_part

/**
 *  \brief v3 -> v1 - v2
 *
 *  Subtracts the input vector @v2 from @v1 and places the result in @v3.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 */
void vec_sub(vector_t *v1, vector_t *v2, vector_t *v3);
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
void vec_sub_part(vector_t *v1, vector_t *v2, vector_t *v3,
                  unsigned long start, unsigned long end);
#define spx_vec_sub_part vec_sub_part

/**
 *  Returns the product of the input vectors @v1 and @v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @return                 the product of the input vectors.
 */
double vec_mul(const vector_t *v1, const vector_t *v2);
#define spx_vec_mul vec_mul

/**
 *  \brief Returns v1[start...end-1] * v2[start...end-1]
 *
 *  Returns the product of the part [@start, @end) of the input vectors
 *  @v1 and @v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 *  @return                 the product of the input vectors.
 */
double vec_mul_part(const vector_t *v1, const vector_t *v2,
                    unsigned long start, unsigned long end);
#define spx_vec_mul_part vec_mul_part

/**
 *  Reorders the input vector @v according to the permutation @p, leaving
 *  the original vector intact.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] p            a permutation.
 *  @return                 the permuted input vector.
 */
spx_error_t spx_vec_reorder(vector_t *v, perm_t *p);

/**
 *  Inverse-reorders the permuted input vector @v1, according
 *  to the permutation @p, storing the result in @v2.
 *
 *  @param[in] v1           a valid permuted vector object.
 *  @param[in/out] v2       the output inverse-reordered vector object.
 *  @param[in] p            a permutation.
 *  @return                 an error code.
 */
spx_error_t spx_vec_inv_reorder(vector_t *v, perm_t *p);

/**
 *  Copies the elements of @v1 to @v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
void vec_copy(const vector_t *v1, vector_t *v2);
#define spx_vec_copy vec_copy

/**
 *  Compares the elements of @v1 and @v2. If they are equal it returns 0,
 *  else -1.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
int vec_compare(const vector_t *v1, const vector_t *v2);
#define spx_vec_compare vec_compare

/**
 *  Prints the input vector @v.
 *
 *  @param[in] v            a valid vector object.
 *  @return                 error code.
 */
void vec_print(const vector_t *v);
#define spx_vec_print vec_print

/**
 *  Destroys the input vector @v.
 *
 *  @param[in] v            a valid vector object.
 */
void vec_destroy(vector_t *v);
#define spx_vec_destroy vec_destroy

#endif // SPARSEX_MAT_VEC_H

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
