/*
 * vector.h
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "macros.h"

/* helpers for cpp abuse */
#ifndef ELEM_TYPE
#define ELEM_TYPE double
#endif

struct map {
    unsigned int length;
    unsigned int *cpus;
    unsigned int *elems_pos;
};
typedef struct map map_t;

#define _VNAME(elem_type, name) vector_ ## elem_type ## _ ## name
#define _VTYPE(elem_type) _VNAME(elem_type, t)

#define DECLARE_VECTOR(vtype) \
typedef struct { \
	vtype *elements; \
	unsigned long size; \
    int alloc_type;     \
} _VTYPE(vtype); \
\
_VTYPE(vtype) *_VNAME(vtype,create)(unsigned long size); \
_VTYPE(vtype) *_VNAME(vtype,create_from_buff)(ELEM_TYPE *buff, unsigned long size); \
_VTYPE(vtype) *_VNAME(vtype,create_onnode)(unsigned long size, int node); \
_VTYPE(vtype) *_VNAME(vtype,create_interleaved)(unsigned long size, size_t *parts, int nr_parts, const int *nodes); \
void _VNAME(vtype,destroy)(_VTYPE(vtype) *v); \
void _VNAME(vtype,init)(_VTYPE(vtype) *v, vtype val); \
void _VNAME(vtype,init_part)(_VTYPE(vtype) *v, vtype val, unsigned long start, unsigned long end); \
void _VNAME(vtype,init_rand_range)(_VTYPE(vtype) *v, vtype max, vtype min); \
void _VNAME(vtype,init_from_map)(_VTYPE(vtype) **v, vtype val, map_t *map); \
void _VNAME(vtype,add)(_VTYPE(vtype) *v1, _VTYPE(vtype) *v2, _VTYPE(vtype) *v3); \
void _VNAME(vtype,add_part)(_VTYPE(vtype) *v1, _VTYPE(vtype) *v2, _VTYPE(vtype) *v3, unsigned long start, unsigned long end); \
void _VNAME(vtype,addmap)(_VTYPE(vtype) *v1, _VTYPE(vtype) **v2, _VTYPE(vtype) *v3, map_t *map); \
int _VNAME(vtype,compare)(_VTYPE(vtype) *v1, _VTYPE(vtype) *v2); \
void _VNAME(vtype,print)(_VTYPE(vtype) *v);

DECLARE_VECTOR(float)
DECLARE_VECTOR(double)

#undef _VNAME
#undef _VTYPE

#define VECTOR_NAME(name) CON3(vector_, ELEM_TYPE, name)
#define VECTOR_TYPE VECTOR_NAME(_t)

#endif /* __VECTOR_H__ */

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
