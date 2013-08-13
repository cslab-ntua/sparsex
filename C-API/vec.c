/**
 * libcsx/vec.c -- Vector interface.
 *                 Essentially wrappers around the vector interface in lib/spm.
 *
 * Copyright (C) 2013, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2013, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "vec.h"

vector_t *vec_create(unsigned long size)
{
    vector_t *v = INVALID_VEC;
    v = (vector_t *) VECTOR_NAME(_create)(size);

    return v;
}

vector_t *vec_create_from_buff(value_t *x, unsigned long size)
{
    /* Check validity of input argument */
    if (!x) {
        SETERROR_1(LIBCSX_ERR_ARG_INVALID, "invalid input array");
        return INVALID_VEC;
    }

    vector_t *v = INVALID_VEC;
	v = (vector_t *) VECTOR_NAME(_create_from_buff)(x, size);

    return v;
}

vector_t *vec_create_random(unsigned long size)
{
    VECTOR_TYPE *x = NULL;

    x = VECTOR_NAME(_create)(size);
    VECTOR_NAME(_init_rand_range)(x, (ELEM_TYPE) -0.001, (ELEM_TYPE) 0.001);

    return (vector_t *) x;
}

libcsx_error_t vec_set_entry(vector_t *v, index_t idx, value_t val)
{
    /* Check if index is out of bounds */
    if (idx <= 0 || idx > v->size) {
        SETERROR_0(LIBCSX_OUT_OF_BOUNDS);
        return LIBCSX_OUT_OF_BOUNDS;
    }

    v->elements[idx - 1] = val;

    return LIBCSX_SUCCESS;
}

vector_t *vec_scale(const vector_t *v, value_t num)
{
    VECTOR_TYPE *scaled_v = NULL;

    scaled_v = VECTOR_NAME(_create)(v->size);
    VECTOR_NAME(_scale)((VECTOR_TYPE *) v, scaled_v, num);

    return (vector_t *) scaled_v;
}

vector_t *vec_scale_original(vector_t *v, value_t num)
{
    VECTOR_NAME(_scale)((VECTOR_TYPE *) v, (VECTOR_TYPE *) v, num);
    return (vector_t *) v;
}

vector_t *vec_add(vector_t *v1, const vector_t *v2)
{
    VECTOR_NAME(_add)((VECTOR_TYPE *) v1, (VECTOR_TYPE *) v2,
                      (VECTOR_TYPE *) v1);

    return (vector_t *) v1;
}

vector_t *vec_reorder(const vector_t *v, perm_t *p)
{
    index_t i;
    vector_t *permuted_v = NULL;

    permuted_v = VECTOR_NAME(_create)(v->size);
    for (i = 0; i < v->size; i++) {
        permuted_v->elements[p[i]] = v->elements[i];
    }

    return permuted_v;
}

vector_t *vec_inv_reorder(vector_t *v, perm_t *p)
{
    index_t i;
    vector_t *permuted_v = NULL;

    permuted_v = VECTOR_NAME(_create)(v->size);
    for (i = 0; i < v->size; i++) {
        permuted_v->elements[i] = v->elements[p[i]];
    }

    return permuted_v;
}

libcsx_error_t vec_print(vector_t *v)
{
    VECTOR_NAME(_print)((VECTOR_TYPE *) v);

    return LIBCSX_SUCCESS;
}

libcsx_error_t vec_destroy(vector_t *v)
{
    VECTOR_NAME(_destroy)((VECTOR_TYPE *) v);

    return LIBCSX_SUCCESS;
}
