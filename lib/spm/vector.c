/*
 * vector.c
 *
 * Copyright (C) 2007-2011, Computing Systems Laboratory (CSLab), NTUA
 * Copyright (C) 2007-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <numa.h>
#include <numaif.h>
#include <sys/mman.h>

#include "vector.h"

enum {
    ALLOC_STD = 1,
    ALLOC_NUMA,
    ALLOC_MMAP,
    ALLOC_OTHER,
};

/*
 * vector: align allocation to 12(=3x4) elements.
 * This (quick-)fixes a bug in bcrs spmv code, when the row block size
 * does not align with the y vector. As a results writes are
 * perfomed beyond the end of the y vector
 */
#define ALIGN_BOUND (12)
#define ALIGN(ul)  (ul + (ul % ALIGN_BOUND))

/*
 *	Create a vector from an input buffer. The elements are not copied.
 */ 
VECTOR_TYPE *VECTOR_NAME(_create_from_buff)(ELEM_TYPE *buff, unsigned long size)
{
	VECTOR_TYPE *v = malloc(sizeof(VECTOR_TYPE));
	if (!v) {
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	v->size = size;
	v->alloc_type = ALLOC_OTHER;
	v->elements = buff;
	return v;
}

VECTOR_TYPE *VECTOR_NAME(_create)(unsigned long size)
{
	VECTOR_TYPE *v = malloc(sizeof(VECTOR_TYPE));
	if ( !v ) {
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	v->size = size;
    v->alloc_type = ALLOC_STD;
	v->elements = malloc(sizeof(ELEM_TYPE)*(size + 12));
	if ( !v->elements ){
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	return v;
}

VECTOR_TYPE *VECTOR_NAME(_create_onnode)(unsigned long size, int node)
{
    VECTOR_TYPE *v = numa_alloc_onnode(sizeof(VECTOR_TYPE), node);
    if (!v) {
        perror("numa_alloc_onnode");
        exit(1);
    }

    v->size = size;
    v->alloc_type = ALLOC_NUMA;
    v->elements = numa_alloc_onnode(sizeof(ELEM_TYPE)*size, node);
    if (!v->elements) {
        perror("numa_alloc_onnode");
        exit(1);
    }

    return v;
}

VECTOR_TYPE *VECTOR_NAME(_create_interleaved)(unsigned long size,
                                              size_t *parts,
                                              int nr_parts,
                                              const int *nodes) {
    VECTOR_TYPE *v = mmap(NULL, sizeof(VECTOR_TYPE),
                          PROT_READ | PROT_WRITE,
                          MAP_ANONYMOUS | MAP_PRIVATE, 0, 0);
    if (v == (void *) -1) {
        perror("mmap");
        exit(1);
    }

    v->size = size;
    v->alloc_type = ALLOC_MMAP;
    v->elements = mmap(NULL, sizeof(ELEM_TYPE)*size,
                       PROT_READ | PROT_WRITE,
                       MAP_ANONYMOUS | MAP_PRIVATE, 0, 0);
    if (v->elements == (void *) -1) {
        perror("mmap");
        exit(1);
    }

    int pagesize = numa_pagesize();
    int nodes_max = numa_max_possible_node();

#define PAGE_ALIGN(addr) (void *)((unsigned long) addr & ~(pagesize-1))
    /*
     * Bind parts to specific nodes
     * All parts must be page aligned
     */
    ELEM_TYPE *curr_part = v->elements;
    int i;
    for (i = 0; i < nr_parts; i++) {
        size_t  part_size = parts[i]*sizeof(ELEM_TYPE);
        size_t  rem = part_size % pagesize;
        while (rem < pagesize / 2 && i < nr_parts - 1) {
            /* Leave the page for the next partition */
            part_size -= sizeof(ELEM_TYPE);
            rem = part_size % pagesize;
        }

        unsigned long nodemask = 1 << nodes[i];
        if (mbind(PAGE_ALIGN(curr_part), part_size,
                  MPOL_BIND, &nodemask, nodes_max, 0) < 0) {
            perror("mbind");
            exit(1);
        }

        parts[i] = part_size / sizeof(ELEM_TYPE);
        curr_part += parts[i];
    }

#undef PAGE_ALIGN
    return v;
}


void VECTOR_NAME(_destroy)(VECTOR_TYPE *v)
{
    if (v->alloc_type == ALLOC_STD) {
        free(v->elements);
        free(v);
    } else if (v->alloc_type == ALLOC_NUMA) {
        numa_free(v->elements, sizeof(ELEM_TYPE)*v->size);
        numa_free(v, sizeof(VECTOR_TYPE));
    } else if (v->alloc_type == ALLOC_MMAP) {
        munmap(v->elements, sizeof(ELEM_TYPE)*v->size);
        munmap(v, sizeof(VECTOR_TYPE));
    } else if (v->alloc_type == ALLOC_OTHER) {
        /* Just free our stuff; elements are supplied from user */
        free(v);
    } else {
        assert(0 && "unknown allocation type");
    }
}

void VECTOR_NAME(_init)(VECTOR_TYPE *v, ELEM_TYPE val)
{
	unsigned long i;
	for (i=0 ; i<v->size; i++){
		v->elements[i] = val;
	}
}

void VECTOR_NAME(_init_rand_range)(VECTOR_TYPE *v, ELEM_TYPE max, ELEM_TYPE min)
{
	unsigned long i;
	ELEM_TYPE val;
	for (i=0; i<v->size; i++){
		val = ((ELEM_TYPE)(rand()+i)/((ELEM_TYPE)RAND_MAX+1));
		v->elements[i] = min + val*(max-min);
	}
}

void VECTOR_NAME(_add)(VECTOR_TYPE *v1, VECTOR_TYPE *v2, VECTOR_TYPE *v3)
{
	unsigned long i;
	
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ", v1->size, v2->size, v3->size);
		exit(1);
	}
	
	for (i = 0; i < v1->size; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

void VECTOR_NAME(_addpart)(VECTOR_TYPE *v1, VECTOR_TYPE *v2, VECTOR_TYPE *v3, unsigned long start, unsigned long end)
{
	unsigned long i;
	
	if (v1->size != v2->size || v1->size != v3->size) {
		fprintf(stderr, "v1->size=%lu v2->size=%lu v3->size=%lu differ", v1->size, v2->size, v3->size);
		exit(1);
	}
	
	if (start > v1->size || end > v1->size || start > end) {
		fprintf(stderr, "start=%lu end=%lu v->size=%lu not compatible", start, end, v1->size);
		exit(1);
	}
	
	for (i = start; i < end; i++)
		v3->elements[i] = v1->elements[i] + v2->elements[i];
}

static inline int elems_neq(ELEM_TYPE a, ELEM_TYPE b)
{
	if ( fabs((double)(a-b)/(double)a)  > 1.e-7 ){
//	if (fabs((double) (a - b)) > 1.e-7) {
		return 1;
	}
	return 0;
}

int VECTOR_NAME(_compare)(VECTOR_TYPE *v1, VECTOR_TYPE *v2)
{
	unsigned long i;
	if ( v1->size != v2->size ){
		fprintf(stderr, "v1->size=%lu v2->size=%lu differ", v1->size, v2->size);
		return -2;
	}
	for (i=0; i<v1->size; i++){
		if ( elems_neq(v1->elements[i], v2->elements[i])){
			fprintf(stderr, "element %ld differs: %10.20lf != %10.20lf\n",
			        i, (double)v1->elements[i], (double)v2->elements[i]);
			return -1;
		}
	}

	return 0;
}

void VECTOR_NAME(_print)(VECTOR_TYPE *v)
{
    unsigned long i;

    printf("[ ");
    for (i = 0; i < v->size; i++)
        printf("%lf ", (double) v->elements[i]);
    printf("]\n");
}

// vim:expandtab:tabstop=8:shiftwidth=4:softtabstop=4
