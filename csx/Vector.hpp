/* -*- C++ -*-
 *
 * Vector.hpp -- Vector interface.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "../api/types.h"
#include "Map.hpp"
#include "numa_util.h"
#include "SpmMt.hpp"

#ifdef __cplusplus

#include <cstdlib>
#include <cassert>
#include <numa.h>
#include <math.h>

extern "C" {
#endif

typedef struct vec {
    value_t *elements;
    unsigned long size;
    int alloc_type;
} vector_t;

typedef index_t perm_t;

vector_t *vec_create(unsigned long size, void *arg);
vector_t *vec_create_from_buff(value_t *buff, unsigned long size, 
                               void *arg);
vector_t *vec_create_onnode(unsigned long size, int node);
vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                 int nr_parts, int *nodes);
vector_t *vec_create_random(unsigned long size, void *A);
void vec_destroy(vector_t *v);
void vec_init(vector_t *v, value_t val);
void init_part(vector_t *v, value_t val, unsigned long start,
               unsigned long end);
void vec_init_from_map(vector_t **v, value_t val, map_t *map);
void vec_init_rand_range(vector_t *v, value_t max, value_t min);
void vec_set_entry(vector_t *v, int idx, value_t val);
void vec_add(vector_t *v1, vector_t *v2, vector_t *v3);
void vec_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                  unsigned long start, unsigned long end);
void vec_add_from_map(vector_t *v1, vector_t **v2, vector_t *v3, 
                      map_t *map);
void vec_sub(vector_t *v1, vector_t *v2, vector_t *v3);
void vec_sub_part(vector_t *v1, vector_t *v2, vector_t *v3,
                  unsigned long start, unsigned long end);
value_t vec_mul(const vector_t *v1, const vector_t *v2);
value_t vec_mul_part(const vector_t *v1, const vector_t *v2,
                    unsigned long start, unsigned long end);
void vec_scale(vector_t *v1, vector_t *v2, scalar_t num);
void vec_scale_part(vector_t *v1, vector_t *v2, scalar_t num,
                    unsigned long start, unsigned long end);
void vec_scale_add(vector_t *v1, vector_t *v2, vector_t *v3,
                   scalar_t num);
void vec_scale_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                        double num, unsigned long start,
                        unsigned long end);
void vec_copy(const vector_t *v1, vector_t *v2);
int vec_compare(const vector_t *v1, const vector_t *v2);
vector_t *vec_reorder(const vector_t *v, perm_t *p);
vector_t *vec_inv_reorder(const vector_t *v, perm_t *p);
void vec_print(const vector_t *v);

#ifdef __cplusplus
}
#endif

#endif // VECTOR_HPP
