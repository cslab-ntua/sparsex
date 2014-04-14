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

#include "Map.hpp"
#include "numa_util.h"
#include "SpmMt.hpp"
#include "../api/types.h"

BEGIN_C_DECLS

typedef struct vec {
    spx_value_t *elements;
    unsigned long size;
    int alloc_type;
} spx_vector_t;

typedef spx_index_t spx_perm_t;

spx_vector_t *vec_create(unsigned long size, void *arg);
spx_vector_t *vec_create_from_buff(spx_value_t *buff, unsigned long size, 
                                   void *arg);
spx_vector_t *vec_create_onnode(unsigned long size, int node);
spx_vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                     int nr_parts, int *nodes);
spx_vector_t *vec_create_random(unsigned long size, void *A);
void vec_destroy(spx_vector_t *v);
void vec_init(spx_vector_t *v, spx_value_t val);
void init_part(spx_vector_t *v, spx_value_t val, spx_index_t start,
               spx_index_t end);
void vec_init_from_map(spx_vector_t **v, spx_value_t val, map_t *map);
void vec_init_rand_range(spx_vector_t *v, spx_value_t max, spx_value_t min);
void vec_set_entry(spx_vector_t *v, int idx, spx_value_t val);
void vec_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3);
void vec_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                  spx_index_t start, spx_index_t end);
void vec_add_from_map(spx_vector_t *v1, spx_vector_t **v2, spx_vector_t *v3, 
                      map_t *map);
void vec_sub(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3);
void vec_sub_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                  spx_index_t start, spx_index_t end);
spx_value_t vec_mul(const spx_vector_t *v1, const spx_vector_t *v2);
spx_value_t vec_mul_part(const spx_vector_t *v1, const spx_vector_t *v2,
                         spx_index_t start, spx_index_t end);
void vec_scale(spx_vector_t *v1, spx_vector_t *v2, spx_scalar_t num);
void vec_scale_part(spx_vector_t *v1, spx_vector_t *v2, spx_scalar_t num,
                    spx_index_t start, spx_index_t end);
void vec_scale_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                   spx_scalar_t num);
void vec_scale_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                        spx_scalar_t num, spx_index_t start, spx_index_t end);
void vec_copy(const spx_vector_t *v1, spx_vector_t *v2);
int vec_compare(const spx_vector_t *v1, const spx_vector_t *v2);
spx_vector_t *vec_reorder(const spx_vector_t *v, spx_perm_t *p);
spx_vector_t *vec_inv_reorder(const spx_vector_t *v, spx_perm_t *p);
void vec_print(const spx_vector_t *v);

END_C_DECLS

#endif // VECTOR_HPP
