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

#ifndef SPARSEX_INTERNALS_VECTOR_HPP
#define SPARSEX_INTERNALS_VECTOR_HPP

#include "sparsex/internals/Map.hpp"
#include "sparsex/internals/numa_util.h"
#include "sparsex/types.h"

BEGIN_C_DECLS

typedef struct vec {
    spx_value_t *elements;
    spx_value_t *ptr_buff;
    size_t size;
    int alloc_type;
    int copy_mode;
} vector_t;

vector_t *VecCreate(size_t size);
vector_t *VecCreateFromBuff(spx_value_t *buff, size_t size, int mode);
vector_t *VecCreateOnnode(size_t size, int node);
vector_t *VecCreateInterleaved(size_t size, size_t *parts, int nr_parts,
                               int *nodes);
vector_t *VecCreateRandom(size_t size);
void VecDestroy(vector_t *v);
void VecInit(vector_t *v, spx_value_t val);
void VecInitPart(vector_t *v, spx_value_t val, spx_index_t start,
                 spx_index_t end);
void VecInitFromMap(vector_t **v, spx_value_t val, map_t *map);
void VecInitRandRange(vector_t *v, spx_value_t max, spx_value_t min);
void VecSetEntry(vector_t *v, spx_index_t idx, spx_value_t val);
void VecAdd(vector_t *v1, vector_t *v2, vector_t *v3);
void VecAddPart(vector_t *v1, vector_t *v2, vector_t *v3, spx_index_t start,
                spx_index_t end);
void VecAddFromMap(vector_t *v1, vector_t **v2, vector_t *v3, map_t *map);
void VecSub(vector_t *v1, vector_t *v2, vector_t *v3);
void VecSubPart(vector_t *v1, vector_t *v2, vector_t *v3, spx_index_t start,
                spx_index_t end);
spx_value_t VecMult(const vector_t *v1, const vector_t *v2);
spx_value_t VecMultPart(const vector_t *v1, const vector_t *v2,
                        spx_index_t start, spx_index_t end);
void VecScale(vector_t *v1, vector_t *v2, spx_scalar_t num);
void VecScalePart(vector_t *v1, vector_t *v2, spx_scalar_t num,
                  spx_index_t start, spx_index_t end);
void VecScaleAdd(vector_t *v1, vector_t *v2, vector_t *v3, spx_scalar_t num);
void VecScaleAddPart(vector_t *v1, vector_t *v2, vector_t *v3, spx_scalar_t num,
                     spx_index_t start, spx_index_t end);
void VecCopy(const vector_t *v1, vector_t *v2);
int VecCompare(const vector_t *v1, const vector_t *v2);
void VecPrint(const vector_t *v);

END_C_DECLS

#endif  // SPARSEX_INTERNALS_VECTOR_HPP
