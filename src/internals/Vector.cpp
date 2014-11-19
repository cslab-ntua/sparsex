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

#include <sparsex/internals/Affinity.hpp>
#include <sparsex/internals/Config.hpp>
#include <sparsex/internals/Vector.hpp>
#include <sparsex/internals/logger/Logger.hpp>
#include <cstdlib>
#include <cstdarg>
#include <cassert>
#include <numa.h>
#include <math.h>
#include <sys/mman.h>
#include <omp.h>

using namespace std;
using namespace sparsex::runtime;

namespace internal {

enum alloc_type {
	ALLOC_STD = 1,
	ALLOC_NUMA,
	ALLOC_MMAP,
	ALLOC_OTHER
};

enum vec_mode {
	VEC_MODE_AS_IS = 43,
	VEC_MODE_TUNED,
    VEC_MODE_INVALID
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
		LOG_ERROR << "malloc failed\n";
		exit(1);
	}

	v->elements = (spx_value_t *) calloc(size, sizeof(spx_value_t));
	if (!v->elements) {
		LOG_ERROR << "malloc failed\n";
		exit(1);
	}

	v->size = size;
	v->alloc_type = internal::ALLOC_STD;
	v->vec_mode = internal::VEC_MODE_INVALID;

	return v;
}

vector_t *VecCreateOnnode(size_t size, int node)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), node);
	v->elements = (spx_value_t *) alloc_onnode(sizeof(spx_value_t)*size, node);
	if (!v->elements) {
		LOG_ERROR << "allocation failed\n";
		exit(1);
	}

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->vec_mode = internal::VEC_MODE_INVALID;

	return v;
}

vector_t *VecCreateInterleaved(size_t size, size_t nr_parts, size_t *parts,
                               int *nodes, spx_index_t *row_start,
                               spx_index_t *row_end)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), nodes[0]);
	v->elements = (spx_value_t *) alloc_interleaved(size*sizeof(*v->elements),
                                                    parts, nr_parts, nodes);
	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->vec_mode = internal::VEC_MODE_INVALID;
    VecInit(v, 0);

	return v;
}

vector_t *VecCreateFromBuff(spx_value_t *buff, size_t size, size_t nr_parts,
                            int *cpus, spx_index_t *row_start,
                            spx_index_t *row_end, int mode)
{
	vector_t *v = (vector_t *) malloc(sizeof(vector_t));
	if (!v) {
		LOG_ERROR << "malloc failed\n";
		exit(1);
	}

    if (mode == internal::VEC_MODE_AS_IS) {
        // We don't know how the data was allocated
        v->vec_mode = internal::VEC_MODE_AS_IS;
        v->elements = buff;
    } else if (mode == internal::VEC_MODE_TUNED) {
        v->vec_mode = internal::VEC_MODE_TUNED;
#if SPX_USE_NUMA
        v->elements = (spx_value_t *) malloc(sizeof(spx_value_t)*size);
        if (!v->elements) {
            LOG_ERROR << "malloc failed\n";
            exit(1);
        }

        // Perform parallel copy
        omp_set_num_threads(nr_parts);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            setaffinity_oncpu(cpus[tid]);
            spx_index_t start = row_start[tid];
            spx_index_t end = row_end[tid];
            for (spx_index_t i = start; i < end; i++)
                v->elements[i] = buff[i];
        }
#else
        v->elements = buff;
#endif
    } else {
        LOG_ERROR << "invalid vector mode";
        exit(1);
    }

    v->alloc_type = internal::ALLOC_OTHER;
	v->size = size;

	return v;
}

vector_t *VecCreateRandom(size_t size)
{
    vector_t *v = VecCreate(size);
    VecInitRandRange(v, (spx_value_t) -0.1, (spx_value_t) 0.1);

    return v;
}

vector_t *VecCreateRandomInterleaved(size_t size, size_t nr_parts,
                                     size_t *parts, int *nodes,
                                     spx_index_t *row_start,
                                     spx_index_t *row_end)
{
    // Save parts since they are modified by alloc_interleaved
    size_t *tmp = new size_t[nr_parts];
    memcpy(tmp, parts, nr_parts * sizeof(size_t));
	vector_t *v = VecCreateInterleaved(size, nr_parts, tmp, nodes,
                                       row_start, row_end);
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

void VecInitRandRangePart(vector_t *v, spx_index_t start, spx_index_t end,
                          spx_value_t max, spx_value_t min)
{
	spx_value_t val;

	for (spx_index_t i = start; i < end; i++) {
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

#if SPX_USE_NUMA
    for (size_t i = 0; i < v1->size; i++)
        v2->elements[i] = v1->elements[i];
#else
    if (v1->vec_mode == internal::VEC_MODE_INVALID &&
        v2->vec_mode == internal::VEC_MODE_INVALID)
        memcpy(v2->elements, v1->elements, v1->size*sizeof(spx_value_t));
    else
        for (size_t i = 0; i < v1->size; i++)
            v2->elements[i] = v1->elements[i];
#endif
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
