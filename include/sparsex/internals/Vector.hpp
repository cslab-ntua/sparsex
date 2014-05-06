/*
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

#include "sparsex/internals/Map.hpp"
#include "sparsex/internals/numa_util.h"
#include "sparsex/internals/SpmMt.hpp"
#include "sparsex/types.h"

BEGIN_C_DECLS

typedef struct vec {
    spx_value_t *elements;
    spx_value_t *ptr_buff;
    unsigned long size;
    int alloc_type;
    int copy_mode;
} vector_t;

vector_t *vec_create(unsigned long size, void *arg);
vector_t *vec_create_from_buff(spx_value_t *buff, unsigned long size, 
                               void *arg, int mode);
vector_t *vec_create_onnode(unsigned long size, int node);
vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                 int nr_parts, int *nodes);
vector_t *vec_create_random(unsigned long size, void *A);
void spx_vec_destroy(vector_t *v);
void spx_vec_init(vector_t *v, spx_value_t val);
void spx_vec_init_part(vector_t *v, spx_value_t val, spx_index_t start,
                       spx_index_t end);
void spx_vec_init_from_map(vector_t **v, spx_value_t val, map_t *map);
void spx_vec_init_rand_range(vector_t *v, spx_value_t max, spx_value_t min);
void vec_set_entry(vector_t *v, int idx, spx_value_t val);
void spx_vec_add(vector_t *v1, vector_t *v2, vector_t *v3);
void spx_vec_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                      spx_index_t start, spx_index_t end);
void spx_vec_add_from_map(vector_t *v1, vector_t **v2, vector_t *v3, 
                          map_t *map);
void spx_vec_sub(vector_t *v1, vector_t *v2, vector_t *v3);
void spx_vec_sub_part(vector_t *v1, vector_t *v2, vector_t *v3,
                      spx_index_t start, spx_index_t end);
spx_value_t spx_vec_mul(const vector_t *v1, const vector_t *v2);
spx_value_t spx_vec_mul_part(const vector_t *v1, const vector_t *v2,
                             spx_index_t start, spx_index_t end);
void spx_vec_scale(vector_t *v1, vector_t *v2, spx_scalar_t num);
void spx_vec_scale_part(vector_t *v1, vector_t *v2, spx_scalar_t num,
                        spx_index_t start, spx_index_t end);
void spx_vec_scale_add(vector_t *v1, vector_t *v2, vector_t *v3,
                       spx_scalar_t num);
void spx_vec_scale_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                            spx_scalar_t num, spx_index_t start,
                            spx_index_t end);
void spx_vec_copy(const vector_t *v1, vector_t *v2);
int spx_vec_compare(const vector_t *v1, const vector_t *v2);
void spx_vec_print(const vector_t *v);

END_C_DECLS

#endif // VECTOR_HPP
