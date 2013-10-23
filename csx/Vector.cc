/* -*- C++ -*-
 *
 * Vector.cc -- Vector interface.
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

vector_t *vec_create(unsigned long size)
{
	vector_t *v = (vector_t *) malloc(sizeof(vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	v->size = size;
	v->alloc_type = internal::ALLOC_STD;
	v->elements = (double *) malloc(sizeof(double)*(size + 12));
	if (!v->elements) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}

	return v;
}

vector_t *vec_create_from_buff(double *buff, unsigned long size)
{
	vector_t *v = (vector_t *) malloc(sizeof(vector_t));
	if (!v) {
		LOG_ERROR << "malloc\n";
		exit(1);
	}
	v->size = size;
	v->alloc_type = internal::ALLOC_OTHER;
	v->elements = buff;

	return v;
}

vector_t *vec_create_onnode(unsigned long size, int node)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), node);

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (double *) alloc_onnode(sizeof(double)*size, node);

	return v;
}

vector_t *vec_create_interleaved(unsigned long size, size_t *parts,
                                 int nr_parts, int *nodes)
{
	vector_t *v = (vector_t *) alloc_onnode(sizeof(vector_t), nodes[0]);

	v->size = size;
	v->alloc_type = internal::ALLOC_MMAP;
	v->elements = (double *) alloc_interleaved(size*sizeof(*v->elements),
                                               parts, nr_parts, nodes);
	return v;
}

vector_t *vec_create_random(unsigned long size)
{
    vector_t *x = NULL;

    x = vec_create(size);
    vec_init_rand_range(x, (double) -0.001, (double) 0.001);

    return x;
}

void vec_destroy(vector_t *v)
{
	if (v->alloc_type == internal::ALLOC_STD) {
		free(v->elements);
		free(v);
	} else if (v->alloc_type == internal::ALLOC_NUMA) {
		numa_free(v->elements, sizeof(double)*v->size);
		numa_free(v, sizeof(vector_t));
	} else if (v->alloc_type == internal::ALLOC_MMAP) {
		munmap(v->elements, sizeof(double)*v->size);
		munmap(v, sizeof(vector_t));
	} else if (v->alloc_type == internal::ALLOC_OTHER) {
		/* Just free our stuff; elements are supplied from user */
		free(v);
	} else {
		assert(0 && "unknown allocation type");
	}
}

void vec_init(vector_t *v, double val)
{
	for (unsigned long i = 0 ; i < v->size; i++)
		v->elements[i] = val;
}

void vec_init_part(vector_t *v, double val, unsigned long start,
                   unsigned long end)
{
	for (unsigned long i = start; i < end; i++)
		v->elements[i] = val;
}

void vec_init_from_map(vector_t **v, double val, map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (unsigned int i = 0; i < map->length; i++)
		v[cpus[i]]->elements[pos[i]] = 0;
}

void vec_init_rand_range(vector_t *v, double max, double min)
{
	double val;

	for (unsigned long i = 0; i < v->size; i++) {
		val = ((double) (rand()+i) / ((double) RAND_MAX + 1));
		v->elements[i] = min + val*(max-min);
	}
}

void vec_set_entry(vector_t *v, int idx, double val)
{
    /* Check if index is out of bounds */
    if (idx <= 0 || static_cast<unsigned long>(idx) > v->size) {
        exit(1);
    }

    v->elements[idx - 1] = val;
}

void vec_add(vector_t *v1, vector_t *v2, vector_t *v3)
{
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ\n",
		        v1->size, v2->size, v3->size);
		exit(1);
	}

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void vec_add_part(vector_t *v1, vector_t *v2, vector_t *v3, unsigned long start,
                  unsigned long end)
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

	for (unsigned long i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void vec_add_from_map(vector_t *v1, vector_t **v2, vector_t *v3, map_t *map)
{
	unsigned int *cpus = map->cpus;
	unsigned int *pos = map->elems_pos;

	for (unsigned int i = 0; i < map->length; i++)
		v3->elements[pos[i]] =
		    v1->elements[pos[i]] + v2[cpus[i]]->elements[pos[i]];
}

void vec_sub(vector_t *v1, vector_t *v2, vector_t *v3)
{
	assert(v1->size == v2->size && v1->size == v3->size && 
	       "vectors for sub have different size");

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

void vec_sub_part(vector_t *v1, vector_t *v2, vector_t *v3, unsigned long start,
                  unsigned long end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for sub have different size");

	for (unsigned long i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] - v2->elements[i];
}

double vec_mul(const vector_t *v1, const vector_t *v2)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

    double ret = 0;
    for (unsigned long i = 0; i < v1->size; i++)
        ret += v1->elements[i] * v2->elements[i];

	return ret;
}

double vec_mul_part(const vector_t *v1, const vector_t *v2,
                    unsigned long start, unsigned long end)
{
	assert(v1->size == v2->size &&  "vectors for mul have incompatible sizes");

	double ret = 0;
	for (unsigned long i = start; i < end; i++)
		ret += v1->elements[i] * v2->elements[i];

	return ret;
}

void vec_scale(vector_t *v1, vector_t *v2, double num)
{
	assert(v1->size == v2->size && "incompatible vector sizes");

	for (unsigned long i = 0; i < v1->size; i++)
		v2->elements[i] = num * v1->elements[i];
}

void vec_scale_part(vector_t *v1, vector_t *v2, double num,
                    unsigned long start, unsigned long end)
{
	assert(v1->size == v2->size &&  "vectors have incompatible sizes");

	for (unsigned long i = start; i < end; i++)
		v2->elements[i] = num * v1->elements[i];
}

void vec_scale_add(vector_t *v1, vector_t *v2, vector_t *v3, double num)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale have incompatible sizes");

	for (unsigned long i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void vec_scale_add_part(vector_t *v1, vector_t *v2, vector_t *v3, double num,
                        unsigned long start, unsigned long end)
{
	assert(v1->size == v2->size && v1->size == v3->size &&
	       "vectors for scale add have incompatible  sizes");

	for (unsigned long i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + num * v2->elements[i];
}

void vec_copy(const vector_t *v1, vector_t *v2)
{
    assert(v1->size == v2->size && "vectors for copy have different size");

    for (unsigned long i = 0; i < v1->size; i++)
        v2->elements[i] = v1->elements[i];
}

static inline int elems_neq(double a, double b)
{
	if (fabs((double) (a - b) / (double) a)  > 1.e-7)
		return 1;
	return 0;
}

int vec_compare(const vector_t *v1, const vector_t *v2)
{
	if (v1->size != v2->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu differ\n", v1->size,
		        v2->size);
		return -2;
	}
	
	for (unsigned long i=0; i<v1->size; i++) {
		if (elems_neq(v1->elements[i], v2->elements[i])) {
			fprintf(stderr, "element %ld differs: %10.20lf != %10.20lf\n", i,
			        (double) v1->elements[i], (double) v2->elements[i]);
			return -1;
		}
	}

	return 0;
}

vector_t *vec_reorder(const vector_t *v, perm_t *p)
{
    unsigned long i;
    vector_t *permuted_v = NULL;

    permuted_v = vec_create(v->size);
    for (i = 0; i < v->size; i++) {
        permuted_v->elements[p[i]] = v->elements[i];
    }

    return permuted_v;
}

vector_t *vec_inv_reorder(const vector_t *v, perm_t *p)
{
    unsigned long i;
    vector_t *permuted_v = NULL;

    permuted_v = vec_create(v->size);
    for (i = 0; i < v->size; i++) {
        permuted_v->elements[i] = v->elements[p[i]];
    }

    return permuted_v;
}

void vec_print(const vector_t *v)
{
    std::cout << "[ ";
	for (unsigned long i = 0; i < v->size; i++)
        std::cout << v->elements[i] << " ";
    std::cout << "]" << std::endl;
}
