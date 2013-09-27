/**
 * libcsx/vec.h -- Vector interface.
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
 *  \brief Dense array object that represents a permutation.
 */
typedef index_t perm_t;

static inline int
check_vec_dim(const vector_t *x, unsigned long dim)
{
    return (x->size == dim);
}

/**
 *  Creates and returns a valid vector object, whose values must be explicitly
 *  initialized.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @return                 a valid vector object.
 */
vector_t *vec_create(unsigned long size);

/**
 *  Creates and returns a valid vector object, whose values are mapped to a
 *  user-defined array.
 *
 *  @param[in] buff         the user-supplied buffer.
 *  @param[in] size         the size of the buffer.
 *  @return                 a valid vector object.
 */
vector_t *vec_create_from_buff(value_t *buff, unsigned long size);
vector_t *vec_create_onnode(unsigned long size, int node);
vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                 int nr_parts, int *nodes);

/**
 *  Creates and returns a valid vector object, whose values are randomly filled.
 *
 *  @param[in] size         the size of the vector to be created.
 *  @return                 a valid vector object.
 */
vector_t *vec_create_random(unsigned long size);

/**
 *  Initializes the valid vector object @v with @val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 */
void vec_init(vector_t *v, double val);

/**
 *  Initializes the [@start, @end) part of the valid vector object @v
 *  with @val.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] val          the value to fill the vector with.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_init_part(vector_t *v, double val, unsigned long start,
                   unsigned long end);
//void vec_init_from_map(vector_t **v, double val, map_t *map);

/**
 *  Initializes the valid vector object @v with random values in the
 *  range [@min, @max].
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] max          maximum value of initializing range.
 *  @param[in] min          minimum value of initializing range.
 */
void vec_init_rand_range(vector_t *v, double max, double min);

/**
 *  Sets the element at index @idx of vector @v to be equal to @val.
 *  @idx is assumed to be one-based.
 *
 *  @param[in] v            a valid vector object.
 *  @param[in] idx          an index inside the vector.
 *  @param[in] val          the value to be set.
 */
void vec_set_entry(vector_t *v, int idx, double val);

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

/**
 *  \brief v3[start...end-1] -> v1[start...end-1] + num * v2[start...end-1]
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @param[in] v3           a valid vector object.
 *  @param[in] start        starting index.
 *  @param[in] end          ending index.
 */
void vec_scale_add_part(vector_t *v1, vector_t *v2, vector_t *v3, double num,
                        unsigned long start, unsigned long end);

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
//void vec_add_from_map(vector_t *v1, vector_t **v2, vector_t *v3, map_t *map);

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

/**
 *  Returns the product of the input vectors @v1 and @v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 *  @return                 the product of the input vectors.
 */
double vec_mul(const vector_t *v1, const vector_t *v2);

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
vector_t *vec_inv_reorder(const vector_t *v, perm_t *p);

/**
 *  Copies the elements of @v1 to @v2.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
void vec_copy(const vector_t *v1, vector_t *v2);

/**
 *  Compares the elements of @v1 and @v2. If they are equal it returns 0,
 *  else -1.
 *
 *  @param[in] v1           a valid vector object.
 *  @param[in] v2           a valid vector object.
 */
int vec_compare(const vector_t *v1, const vector_t *v2);

/**
 *  Prints the input vector @v.
 *
 *  @param[in] v            a valid vector object.
 *  @return                 error code.
 */
void vec_print(const vector_t *v);

/**
 *  Destroys the input vector @v.
 *
 *  @param[in] v            a valid vector object.
 */
void vec_destroy(vector_t *v);

#endif // LIBCSX_VECTOR_H__
