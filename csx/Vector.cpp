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

#include "Vector.hpp"
#include "Logger.hpp"

#include <sys/mman.h>

namespace internal {
enum alloc_type {
	ALLOC_STD = 1,
	ALLOC_NUMA,
	ALLOC_MMAP,
	ALLOC_OTHER,
};
}

spx_vector_t *vec_create(unsigned long size, void *arg)
{
	spx_vector_t *v = (spx_vector_t *) malloc(sizeof(spx_vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	v->size = size;
	v->alloc_type = internal::ALLOC_STD;
	v->elements = (spx_value_t *) malloc(sizeof(spx_value_t)*(size + 12));
	if (!v->elements) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	return v;
}

spx_vector_t *vec_create_from_buff(spx_value_t *buff, unsigned long size,
                                   void *arg)
{
	spx_vector_t *v = (spx_vector_t *) malloc(sizeof(spx_vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}
	v->size = size;
	v->alloc_type = internal::ALLOC_OTHER;
	v->elements = buff;

	return v;
}

spx_vector_t *vec_create_onnode(unsigned long size, int node)
{
	spx_vector_t *v = (spx_vector_t *) alloc_onnode(sizeof(spx_vector_t), node);

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (spx_value_t *) alloc_onnode(sizeof(spx_value_t)*size, node);

	return v;
}

spx_vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                     int nr_parts, int *nodes)
{
	spx_vector_t *v = (spx_vector_t *) alloc_onnode(sizeof(spx_vector_t),
                                                    nodes[0]);

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (spx_value_t *) alloc_interleaved(size*sizeof(*v->elements),
                                                    parts, nr_parts, nodes);
	return v;
}

spx_vector_t *vec_create_random(unsigned long size, void *A)
{
    spx_vector_t *x = NULL;

    x = vec_create(size, NULL);
    vec_init_rand_range(x, (spx_value_t) -0.1, (spx_value_t) 0.1);

    return x;
}

void vec_destroy(spx_vector_t *v)
{
	if (v->alloc_type == internal::ALLOC_STD) {
		free(v->elements);
		free(v);
	} else if (v->alloc_type == internal::ALLOC_NUMA) {
		numa_free(v->elements, sizeof(spx_value_t)*v->size);
		numa_free(v, sizeof(spx_vector_t));
	} else if (v->alloc_type == internal::ALLOC_MMAP) {
		munmap(v->elements, sizeof(spx_value_t)*v->size);
		munmap(v, sizeof(spx_vector_t));
	} else if (v->alloc_type == internal::ALLOC_OTHER) {
		/* Just free our stuff; elements are supplied from user */
		free(v);
	} else {
		assert(0 && "unknown allocation type");
	}
}

void vec_init(spx_vector_t *v, spx_value_t val)
{
	for (unsigned long i = 0 ; i < v->size; i++)
		v->elements[i] = val;
}

void vec_init_part(spx_vector_t *v, spx_value_t val, spx_index_t start,
                   spx_index_t end)
{
	for (spx_index_t i = start; i < end; i++)
		v->elements[i] = val;
}

void vec_init_from_map(spx_vector_t **v, spx_value_t val, map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (unsigned int i = 0; i < map->length; i++)
		v[cpus[i]]->elements[pos[i]] = val;
}

void vec_init_rand_range(spx_vector_t *v, spx_value_t max, spx_value_t min)
{
	spx_value_t val;

	for (unsigned long i = 0; i < v->size; i++) {
		val = ((spx_value_t) (rand()+i) / ((spx_value_t) RAND_MAX + 1));
		v->elements[i] = min + val*(max-min);
	}
}

void vec_set_entry(spx_vector_t *v, spx_index_t idx, spx_value_t val)
{
    /* Check if index is out of bounds */
    if (idx <= 0 || static_cast<unsigned long>(idx) > v->size) {
        exit(1);
    }

    v->elements[idx - 1] = val;
}

void vec_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void vec_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3, 
                  spx_index_t start, spx_index_t end)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	if (start > v1->size || end > v1->size || start > end) {
		fprintf(stderr, "start=%lu end=%lu v->size=%lu not compatible\n",
		        start, end, v1->size);
		exit(1);
	}

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void vec_add_from_map(spx_vector_t *v1, spx_vector_t **v2, spx_vector_t *v3,
                      map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (unsigned int i = 0; i < map->length; i++)
		v3->elements[pos[i]] =
		    v1->elements[pos[i]] + v2[cpus[i]]->elements[pos[i]];
}

void vec_sub(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3)
{
	assert(v1->size == v2->size && v1->size == v3->size && 
	       "vectors for sub have different size");

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

void vec_sub_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                  spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for sub have different size");

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

spx_value_t vec_mul(const spx_vector_t *v1, const spx_vector_t *v2)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

    spx_value_t ret = 0;
    for (unsigned long i = 0; i < v1->size; i++)
        ret += v1->elements[i] * v2->elements[i];

	return ret;
}

spx_value_t vec_mul_part(const spx_vector_t *v1, const spx_vector_t *v2,
                         spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

	spx_value_t ret = 0;
	for (spx_index_t i = start; i < end; i++)
		ret += v1->elements[i] * v2->elements[i];

	return ret;
}

void vec_scale(spx_vector_t *v1, spx_vector_t *v2, spx_scalar_t num)
{
	assert(v1->size == v2->size && "incompatible vector sizes");

	for (unsigned long i = 0; i < v1->size; i++)
		v2->elements[i] = num * v1->elements[i];
}

void vec_scale_part(spx_vector_t *v1, spx_vector_t *v2, spx_scalar_t num,
                    spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size &&  "vectors have incompatible sizes");

	for (spx_index_t i = start; i < end; i++)
		v2->elements[i] = num * v1->elements[i];
}

void vec_scale_add(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                   spx_scalar_t num)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale have incompatible sizes");

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void vec_scale_add_part(spx_vector_t *v1, spx_vector_t *v2, spx_vector_t *v3,
                        spx_scalar_t num, spx_index_t start, spx_index_t end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale add have incompatible  sizes");

	for (spx_index_t i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void vec_copy(const spx_vector_t *v1, spx_vector_t *v2)
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

int vec_compare(const spx_vector_t *v1, const spx_vector_t *v2)
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

void vec_print(const spx_vector_t *v)
{
    std::cout << "[ ";
	for (unsigned long i = 0; i < v->size; i++)
        std::cout << v->elements[i] << " ";
    std::cout << "]" << std::endl;
}
