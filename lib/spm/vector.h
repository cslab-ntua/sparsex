/*
 * vector.h
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * Copyright (C) 2011,      Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include "macros.h"
#include "numa_util.h"

/* helpers for cpp abuse */
#ifndef ELEM_TYPE
#define ELEM_TYPE double
#endif

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
_VTYPE(vtype) *_VNAME(vtype,create_onnode)(unsigned long size, int node); \
_VTYPE(vtype) *_VNAME(vtype,create_interleaved)(unsigned long size, size_t *parts, int nr_parts, const int *nodes); \
void _VNAME(vtype,destroy)(_VTYPE(vtype) *v); \
void _VNAME(vtype,init)(_VTYPE(vtype) *v, vtype val); \
void _VNAME(vtype,init_part)(_VTYPE(vtype) *v, unsigned long start, unsigned long end, vtype val); \
void _VNAME(vtype,init_rand_range)(_VTYPE(vtype) *v, vtype max, vtype min); \
int _VNAME(vtype,compare)(_VTYPE(vtype) *v1, _VTYPE(vtype) *v2);

DECLARE_VECTOR(float)
DECLARE_VECTOR(double)

#undef _VNAME
#undef _VTYPE

#define VECTOR_NAME(name) CON3(vector_, ELEM_TYPE, name)
#define VECTOR_TYPE VECTOR_NAME(_t)

#endif /* __VECTOR_H__ */
