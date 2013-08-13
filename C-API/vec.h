/**
 * libcsx/vec.h -- Vector interface.
 *                 Essentially wrappers around the vector interface in lib/spm.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef LIBCSX_VECTOR_H__
#define LIBCSX_VECTOR_H__

#include "common.h"
#include "vector.h"

/**
 *  \brief Dense vector object.
 */
typedef VECTOR_TYPE vector_t;

/**
 *  \brief Dense array object that represents a permutation.
 */
typedef index_t perm_t;

static inline int
check_vec_dim(const vector_t *x, index_t dim)
{
    return (x->size == dim);
}

/**
 *  Creates and returns a valid vector object, whose values must be explicitly
 *  filled with the use of #vec_set_entry() or #vec_set_entries().
 *
 *  @param[in] size         the size of the vector to be created.
 *  @return                 a valid vector object.
 */
vector_t *vec_create(unsigned long size);

/**
 *  Creates and returns a valid vector object, whose values are filled by a
 *  user-defined array.
 *
 *  @param[in] buff         the user-supplied buffer.
 *  @param[in] size         the size of the buffer.
 *  @return                 a valid vector object.
 */
vector_t *vec_create_from_buff(value_t *buff, unsigned long size);

/**
 *  Creates and returns a valid vector object, whose values are randomly filled.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @return                 a valid vector object.
 */
vector_t *vec_create_random(unsigned long size);

/**
 *  Sets the element at index @idx of vector @v to be equal to @val.
 *  @idx is assumed to be one-based.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] idx          an index inside the vector.
 *  @param[in] val          the value to be set.
 *  @return                 error code.
 */
libcsx_error_t vec_set_entry(vector_t *v, index_t idx, value_t val);

/**
 *  Scales the input vector @v by a constant value @num and returns
 *  a newly-allocated vector object, leaving the original vector intact.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] num          the const by which to scale the vector.
 *  @return                 a newly-allocated valid vector object.
 */
vector_t *vec_scale(const vector_t *v, value_t num);

/**
 *  Scales the input vector @v by a constant value @num by replacing
 *  the original vector object.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] num          the const by which to scale the vector.
 *  @return                 the input vector object, scaled by "num".
 */
vector_t *vec_scale_original(vector_t *v, value_t num);

/**
 *  Adds the input vectors @v1 and @v2 and places the result in @v1.
 *
 *  @param[in/out] v1       a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @return                 the sum of "v1" and "v2" stored in "v1".
 */
vector_t *vec_add(vector_t *v1, const vector_t *v2);

/**
 *  Reorders the input vector @v according to the permutation @p, leaving
 *  the original vector intact.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] p            a permutation.
 *  @return                 the permuted input vector.
 */
vector_t *vec_reorder(const vector_t *v, perm_t *p);

/**
 *  Restores the permuted input vector @v to its original ordering, according
 *  to the permutation @p, leaving the original vector intact.
 *
 *  @param[in] v            a valid permuted vector object.
 *  @param[in] p            a permutation.
 *  @return                 the permuted input vector in its original ordering.
 */
vector_t *vec_inv_reorder(vector_t *v, perm_t *p);

/**
 *  Prints the input vector @v.
 *
 *  @param[in] v            a valid vector object.
 *  @return                 error code.
 */
libcsx_error_t vec_print(vector_t *v);

/**
 *  Destroys the input vector @v.
 *
 *  @param[in] v            a valid vector object.
 *  @return                 error code.
 */
libcsx_error_t vec_destroy(vector_t *v);

#endif // LIBCSX_VECTOR_H__
