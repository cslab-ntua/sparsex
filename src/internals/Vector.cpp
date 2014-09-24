/*
 * Copyright (C) 2009-2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * Copyright (C) 2012-2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file Vector.cpp
 * \brief Vector interface
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/logger/Logger.hpp>
#include <cstdlib>
#include <cstdarg>
#include <cassert>
#include <numa.h>
#include <math.h>
#include <sys/mman.h>

using namespace std;

namespace internal {

enum alloc_type {
	ALLOC_STD = 1,
	ALLOC_NUMA,
	ALLOC_MMAP,
	ALLOC_OTHER
};

enum copy_mode {
	SHARE = 43,
	COPY
};

template<typename ValueType>
int elems_neq(ValueType a, ValueType b)
{
	if (fabs((ValueType) (a - b) / (ValueType) a)  > 1.e-6)
		return 1;
	return 0;
}

}

vector_t *VecCreate(size_t size)
{
	vector_t *v = (vector_t *) malloc(sizeof(vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	v->size = size;
	v->alloc_type = internal::ALLOC_STD;
	v->elements = (spx_value_t *) calloc(size, sizeof(spx_value_t));
	if (!v->elements) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	return v;
}

vector_t *VecCreateOnnode(size_t size, int node)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), node);
	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (spx_value_t *) alloc_onnode(sizeof(spx_value_t)*size, node);
	if (!v->elements) {
		LOG_ERROR << "allocation failed\n";
		exit(1);
	}

    VecInit(v, 0);
	return v;
}

vector_t *VecCreateInterleaved(size_t size, size_t *parts, int nr_parts,
                               int *nodes)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), nodes[0]);
	v->size = size;
	v->copy_mode = -1;
	v->ptr_buff = NULL;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (spx_value_t *) alloc_interleaved(size*sizeof(*v->elements),
                                                    parts, nr_parts, nodes);
    VecInit(v, 0);
	return v;
}

vector_t *VecCreateFromBuff(spx_value_t *buff, size_t size, int mode)
{
	vector_t *v = (vector_t *) malloc(sizeof(vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}
	v->size = size;
	v->alloc_type = internal::ALLOC_OTHER;
    if (mode == internal::SHARE) {
        v->copy_mode = internal::SHARE;
        v->elements = buff;
    } else if (mode == internal::COPY) {
        v->copy_mode = internal::COPY;
        v->elements = (spx_value_t *) malloc(sizeof(spx_value_t) * size);
        if (!v->elements) {
            LOG_ERROR << "malloc\n";
            exit(1);
        }
        
        memcpy(v->elements, buff, size * sizeof(spx_value_t));
    } else {
        LOG_ERROR << "invalid copy-mode type";
        exit(1);
    }

	return v;
}

vector_t *VecCreateFromBuffInterleaved(spx_value_t *buff, size_t size,
                                       size_t *parts, int nr_parts,
                                       int *nodes, int mode)
{
    // Save parts since they are modified by alloc_interleaved
    size_t *tmp = new size_t[nr_parts];
    memcpy(tmp, parts, nr_parts * sizeof(size_t));
	vector_t *v = VecCreateInterleaved(size, tmp, nr_parts, nodes);
    for (size_t i = 0; i < size; i++)
        v->elements[i] = buff[i];
    v->copy_mode = mode;
    if (mode == internal::SHARE)
        v->ptr_buff = buff;
    else if (mode == internal::COPY)
        v->ptr_buff = NULL;

    print_alloc_status("vector", check_interleaved(v->elements, parts,
                                                   nr_parts,
                                                   nodes));
    delete[] tmp;
    return v;
}

vector_t *VecCreateFromBuffOnnode(spx_value_t *buff, size_t size, int node,
                                  int mode)
{
	vector_t *v = VecCreateOnnode(size, node);
    for (size_t i = 0; i < size; i++)
        v->elements[i] = buff[i];
    v->copy_mode = mode;
    if (mode == internal::SHARE)
        v->ptr_buff = buff;
    else if (mode == internal::COPY)
        v->ptr_buff = NULL;

    return v;
}

vector_t *VecCreateRandom(size_t size)
{
    vector_t *v = VecCreate(size);
    VecInitRandRange(v, (spx_value_t) -0.1, (spx_value_t) 0.1);

    return v;
}

vector_t *VecCreateRandomInterleaved(size_t size, size_t *parts, int nr_parts,
                                     int *nodes)
{
    // Save parts since they are modified by alloc_interleaved
    size_t *tmp = new size_t[nr_parts];
    memcpy(tmp, parts, nr_parts * sizeof(size_t));
	vector_t *v = VecCreateInterleaved(size, tmp, nr_parts, nodes);
    VecInitRandRange(v, (spx_value_t) -0.1, (spx_value_t) 0.1);

    print_alloc_status("vector", check_interleaved(v->elements, tmp,
                                                   nr_parts,
                                                   nodes));
    delete[] tmp;
    return v;
}

void VecDestroy(vector_t *v)
{
	if (v->alloc_type == internal::ALLOC_STD) {
		free(v->elements);
		free(v);
	} else if (v->alloc_type == internal::ALLOC_NUMA) {
		numa_free(v->elements, sizeof(spx_value_t) * v->size);
		numa_free(v, sizeof(vector_t));
	} else if (v->alloc_type == internal::ALLOC_MMAP) {
        if (v->copy_mode == internal::SHARE)
            for (size_t i = 0; i < v->size; i++)
                v->ptr_buff[i] = v->elements[i];
        v->ptr_buff = NULL;
		munmap(v->elements, sizeof(spx_value_t) * v->size);
		munmap(v, sizeof(vector_t));
	} else if (v->alloc_type == internal::ALLOC_OTHER) {
		/* Just free our stuff; elements are supplied from user */
		free(v);
	} else {
		assert(0 && "unknown allocation type");
	}
}

void VecInit(vector_t *v, spx_value_t val)
{
	for (size_t i = 0 ; i < v->size; i++)
		v->elements[i] = val;
}

void VecInitPart(vector_t *v, spx_value_t val, spx_index_t start,
                 spx_index_t end)
{
	for (spx_index_t i = start; i < end; i++)
		v->elements[i] = val;
}

void VecInitFromMap(vector_t **v, spx_value_t val, map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (size_t i = 0; i < map->length; i++)
		v[cpus[i]]->elements[pos[i]] = val;
}

void VecInitRandRange(vector_t *v, spx_value_t max, spx_value_t min)
{
	spx_value_t val;

	for (size_t i = 0; i < v->size; i++) {
		val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
		v->elements[i] = min + val*(max-min);
	}
}

void VecSetEntry(vector_t *v, spx_index_t idx, spx_value_t val)
{
    /* Check if index is out of bounds */
    if (idx <= 0 || static_cast<size_t>(idx) > v->size) {
        exit(1);
    }

    v->elements[idx - 1] = val;
}

void VecAdd(vector_t *v1, vector_t *v2, vector_t *v3)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	for (size_t i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void VecAddPart(vector_t *v1, vector_t *v2, vector_t *v3, spx_index_t start,
                spx_index_t end)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	if (static_cast<size_t>(start) > v1->size ||
        static_cast<size_t>(end) > v1->size || start > end) {
		fprintf(stderr, "start=%lu end=%lu v->size=%lu not compatible\n",
		        static_cast<size_t>(start), static_cast<size_t>(end), v1->size);
		exit(1);
	}

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void VecAddFromMap(vector_t *v1, vector_t **v2, vector_t *v3, map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (size_t i = 0; i < map->length; i++)
		v3->elements[pos[i]] =
		    v1->elements[pos[i]] + v2[cpus[i]]->elements[pos[i]];
}

void VecSub(vector_t *v1, vector_t *v2, vector_t *v3)
{
	assert(v1->size == v2->size && v1->size == v3->size && 
	       "vectors for sub have different size");

	for (size_t i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

void VecSubPart(vector_t *v1, vector_t *v2, vector_t *v3, spx_index_t start,
                spx_index_t end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for sub have different size");

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

spx_value_t VecMult(const vector_t *v1, const vector_t *v2)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

    spx_value_t ret = 0;
    for (size_t i = 0; i < v1->size; i++)
        ret += v1->elements[i] * v2->elements[i];

	return ret;
}

spx_value_t VecMultPart(const vector_t *v1, const vector_t *v2,
                        spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

	spx_value_t ret = 0;
	for (spx_index_t i = start; i < end; i++)
		ret += v1->elements[i] * v2->elements[i];

	return ret;
}

void VecScale(vector_t *v1, vector_t *v2, spx_value_t num)
{
	assert(v1->size == v2->size && "incompatible vector sizes");

	for (size_t i = 0; i < v1->size; i++)
		v2->elements[i] = num * v1->elements[i];
}

void VecScalePart(vector_t *v1, vector_t *v2, spx_value_t num,
                  spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size &&  "vectors have incompatible sizes");

	for (spx_index_t i = start; i < end; i++)
		v2->elements[i] = num * v1->elements[i];
}

void VecScaleAdd(vector_t *v1, vector_t *v2, vector_t *v3, spx_value_t num)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale have incompatible sizes");

	for (size_t i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void VecScaleAddPart(vector_t *v1, vector_t *v2, vector_t *v3,
                     spx_value_t num, spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale add have incompatible  sizes");

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void VecCopy(const vector_t *v1, vector_t *v2)
{
    assert(v1->size == v2->size && "vectors for copy have different size");

    for (size_t i = 0; i < v1->size; i++)
        v2->elements[i] = v1->elements[i];//memcpy
}

int VecCompare(const vector_t *v1, const vector_t *v2)
{
	if (v1->size != v2->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu differ\n", v1->size,
		        v2->size);
		return -2;
	}
	
	for (size_t i = 0; i < v1->size; i++) {
		if (internal::elems_neq(v1->elements[i], v2->elements[i])) {
			fprintf(stderr, "element %ld differs: %10.20f != %10.20f\n", i,
			        v1->elements[i], v2->elements[i]);
			return -1;
		}
	}

	return 0;
}

void VecPrint(const vector_t *v)
{
    cout << "[ ";
	for (size_t i = 0; i < v->size; i++)
        cout << v->elements[i] << " ";
    cout << "]" << endl;
}
