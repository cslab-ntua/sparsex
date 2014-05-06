/*
 * Vector.cpp -- Vector interface.
 *
 * Copyright (C) 2009-2012, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * Copyright (C) 2009-2012, Vasileios Karakasis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

#include "sparsex/internals/Vector.hpp"
#include "sparsex/internals/logger/Logger.hpp"

#include <cstdlib>
#include <cstdarg>
#include <cassert>
#include <numa.h>
#include <math.h>
#include <sys/mman.h>

namespace internal {

enum alloc_type {
	ALLOC_STD = 1,
	ALLOC_NUMA,
	ALLOC_MMAP,
	ALLOC_OTHER,
};

enum copy_mode {
	SHARE = 43,
	COPY
};
}

vector_t *vec_create(unsigned long size, void *arg)
{
	vector_t *v = (vector_t *) malloc(sizeof(vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	v->size = size;
	v->alloc_type = internal::ALLOC_STD;
	v->elements = (spx_value_t *) malloc(sizeof(spx_value_t) * size);
	if (!v->elements) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	return v;
}

vector_t *vec_create_from_buff(spx_value_t *buff, unsigned long size,
                               void *arg, int mode)
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
        LOG_ERROR << "Invalid copy-mode type";
        exit(1);
    }

	return v;
}

vector_t *vec_create_onnode(unsigned long size, int node)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), node);

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (spx_value_t *) alloc_onnode(sizeof(spx_value_t)*size, node);

	return v;
}

vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                 int nr_parts, int *nodes)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t),
                                            nodes[0]);

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (spx_value_t *) alloc_interleaved(size*sizeof(*v->elements),
                                                    parts, nr_parts, nodes);
	return v;
}

vector_t *vec_create_random(unsigned long size, void *A)
{
    vector_t *x = NULL;

    x = vec_create(size, NULL);
    spx_vec_init_rand_range(x, (spx_value_t) -0.1, (spx_value_t) 0.1);

    return x;
}

void spx_vec_destroy(vector_t *v)
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

void spx_vec_init(vector_t *v, spx_value_t val)
{
	for (unsigned long i = 0 ; i < v->size; i++)
		v->elements[i] = val;
}

void spx_vec_init_part(vector_t *v, spx_value_t val, spx_index_t start,
                       spx_index_t end)
{
	for (spx_index_t i = start; i < end; i++)
		v->elements[i] = val;
}

void spx_vec_init_from_map(vector_t **v, spx_value_t val, map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (unsigned int i = 0; i < map->length; i++)
		v[cpus[i]]->elements[pos[i]] = val;
}

void spx_vec_init_rand_range(vector_t *v, spx_value_t max, spx_value_t min)
{
	spx_value_t val;

	for (unsigned long i = 0; i < v->size; i++) {
		val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
		v->elements[i] = min + val*(max-min);
	}
}

void vec_set_entry(vector_t *v, spx_index_t idx, spx_value_t val)
{
    /* Check if index is out of bounds */
    if (idx <= 0 || static_cast<unsigned long>(idx) > v->size) {
        exit(1);
    }

    v->elements[idx - 1] = val;
}

void spx_vec_add(vector_t *v1, vector_t *v2, vector_t *v3)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void spx_vec_add_part(vector_t *v1, vector_t *v2, vector_t *v3, 
                      spx_index_t start, spx_index_t end)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	if ((size_t) start > v1->size || (size_t) end > v1->size || start > end) {
		fprintf(stderr, "start=%lu end=%lu v->size=%lu not compatible\n",
		        (size_t) start, (size_t) end, v1->size);
		exit(1);
	}

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void spx_vec_add_from_map(vector_t *v1, vector_t **v2, vector_t *v3,
                          map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (unsigned int i = 0; i < map->length; i++)
		v3->elements[pos[i]] =
		    v1->elements[pos[i]] + v2[cpus[i]]->elements[pos[i]];
}

void spx_vec_sub(vector_t *v1, vector_t *v2, vector_t *v3)
{
	assert(v1->size == v2->size && v1->size == v3->size && 
	       "vectors for sub have different size");

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

void spx_vec_sub_part(vector_t *v1, vector_t *v2, vector_t *v3,
                      spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for sub have different size");

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

spx_value_t spx_vec_mul(const vector_t *v1, const vector_t *v2)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

    spx_value_t ret = 0;
    for (unsigned long i = 0; i < v1->size; i++)
        ret += v1->elements[i] * v2->elements[i];

	return ret;
}

spx_value_t spx_vec_mul_part(const vector_t *v1, const vector_t *v2,
                             spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

	spx_value_t ret = 0;
	for (spx_index_t i = start; i < end; i++)
		ret += v1->elements[i] * v2->elements[i];

	return ret;
}

void spx_vec_scale(vector_t *v1, vector_t *v2, spx_scalar_t num)
{
	assert(v1->size == v2->size && "incompatible vector sizes");

	for (unsigned long i = 0; i < v1->size; i++)
		v2->elements[i] = num * v1->elements[i];
}

void spx_vec_scale_part(vector_t *v1, vector_t *v2, spx_scalar_t num,
                        spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size &&  "vectors have incompatible sizes");

	for (spx_index_t i = start; i < end; i++)
		v2->elements[i] = num * v1->elements[i];
}

void spx_vec_scale_add(vector_t *v1, vector_t *v2, vector_t *v3,
                       spx_scalar_t num)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale have incompatible sizes");

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void spx_vec_scale_add_part(vector_t *v1, vector_t *v2, vector_t *v3,
                            spx_scalar_t num, spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale add have incompatible  sizes");

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void spx_vec_copy(const vector_t *v1, vector_t *v2)
{
    assert(v1->size == v2->size && "vectors for copy have different size");

    for (unsigned long i = 0; i < v1->size; i++)
        v2->elements[i] = v1->elements[i];
}

static inline int elems_neq(spx_value_t a, spx_value_t b)
{
	if (fabs((spx_value_t) (a - b) / (spx_value_t) a)  > 1.e-7)
		return 1;
	return 0;
}

int spx_vec_compare(const vector_t *v1, const vector_t *v2)
{
	if (v1->size != v2->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu differ\n", v1->size,
		        v2->size);
		return -2;
	}
	
	for (unsigned long i=0; i<v1->size; i++) {
		if (elems_neq(v1->elements[i], v2->elements[i])) {
			fprintf(stderr, "element %ld differs: %10.20lf != %10.20lf\n", i,
			        (spx_value_t) v1->elements[i], (spx_value_t) v2->elements[i]);
			return -1;
		}
	}

	return 0;
}

void spx_vec_print(const vector_t *v)
{
    std::cout << "[ ";
	for (unsigned long i = 0; i < v->size; i++)
        std::cout << v->elements[i] << " ";
    std::cout << "]" << std::endl;
}
